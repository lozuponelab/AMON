import asyncio

import aiohttp

from KEGG_parser.parsers import parse_ko


def get_from_kegg_flat_file(file_loc, list_of_ids=None, parser=parse_ko):
    record_list = list()
    for entry in open(file_loc).read().split('///')[:-1]:
        record = parser(entry)
        if list_of_ids is not None:
            if record['ENTRY'] in list_of_ids:
                record_list.append(record)
        else:
            record_list.append(record)
    return record_list


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


async def download_coroutine(session, url, attempts=10, wait=30):
    for _ in range(attempts):
        async with session.get(url) as response:
            if response.status == 200:
                return await response.text()
            elif response.status != 403:
                raise ValueError('Bad HTTP request status %s: %s\n%s' % (response.status, response.reason, url))
        await asyncio.sleep(wait)

    raise ValueError('KEGG has forbidden request after %s attempts' % attempts)


async def kegg_download_manager(loop, list_of_ids):
    urls = ['http://rest.kegg.jp/get/%s' % '+'.join(chunk) for chunk in chunks(list(list_of_ids), 10)]

    async with aiohttp.ClientSession(loop=loop) as session:
        tasks = [download_coroutine(session, url) for url in urls]
        results = await asyncio.gather(*tasks)

    return [raw_record for raw_records in results for raw_record in raw_records.split('///')[:-1]]


def get_from_kegg_api(loop, list_of_ids, parser):
    return [parser(raw_record) for raw_record in loop.run_until_complete(kegg_download_manager(loop, list_of_ids))]


def get_kegg_record_dict(list_of_ids, parser, records_file_loc=None, verbose=False):
    if records_file_loc is None:
        loop = asyncio.get_event_loop()
        records = get_from_kegg_api(loop, list_of_ids, parser)
    else:
        records = get_from_kegg_flat_file(records_file_loc, list_of_ids, parser)
    if verbose:
        print("%s records acquired" % len(records))
    return {record['ENTRY']: record for record in records}


async def kegg_link_download_manager(loop, link1, link2):
    async with aiohttp.ClientSession(loop=loop) as session:
        url = 'http://rest.kegg.jp/link/%s/%s' % (link1, link2)
        tasks = [download_coroutine(session, url)]
        results = await asyncio.gather(*tasks)
    link_dict = dict()
    if len(results) != 1:
        raise ValueError('Result had more than one value')
    for link in results[0].strip().split('\n'):
        obj_1, obj_2 = link.strip().split()
        obj_1 = obj_1.strip().split(':')[1]
        obj_2 = obj_2.strip().split(':')[1]
        link_dict[obj_2] = obj_1
    return link_dict


def get_kegg_link_from_api(link1, link2):
    loop = asyncio.get_event_loop()
    return loop.run_until_complete(kegg_link_download_manager(loop, link1, link2))
