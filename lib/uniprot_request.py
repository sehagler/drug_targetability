# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:41:23 2023

@author: haglers
"""

#
import json
import re
import requests
from requests.adapters import HTTPAdapter, Retry
import time
from urllib.parse import urlparse, parse_qs, urlencode
from xml.etree import ElementTree
import zlib

#
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

#
retries = \
    Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

#
def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])

#
def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise

#
def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results

#
def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text

#
def divide_list_into_chunks(lst, chunk_length): 
    chunks = []
    if chunk_length is not None:
        stop_flg = False
        while not stop_flg:
            if len(lst) > chunk_length:
                chunks.append(lst[:chunk_length])
                lst = lst[chunk_length:]
            elif len(lst) > 0:
                chunks.append(lst)
                stop_flg = True
            else:
                stop_flg = True
    else:
        chunks.append(lst)
    return chunks

#
def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)

#
def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]

#
def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results

#
def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)

#
def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

#
def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""

#
def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)

#
def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")
    
#
def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(f"{API_URL}/idmapping/run",
                            data={"from": from_db, "to": to_db, "ids": ",".join(ids)})
    check_response(request)
    return request.json()["jobId"]

#
class Uniprot_request_object(object):
    
    #
    def __init__(self):
        self.chunk_length = None
    
    #
    def _map_uniprot_accession_to_gene_name(self, uniprot_accession_list):
        try:
            job_id = submit_id_mapping(
                from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=uniprot_accession_list
            )
            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                results = get_id_mapping_results_search(link)
            request_results_list = []
            for item0 in results['results']:
                for item1 in item0['to']['genes']:
                    request_results_list.append(item1['geneName']['value'])
        except Exception as error:
            print('Uniprot request failed')
            print(error)
            request_results_list = None
        return request_results_list
    
    #
    def map_gene_name_to_uniprot_accession(self, gene_name):
        if isinstance(gene_name, list):
            gene_name_list = gene_name
        else:
            gene_name_list = [ gene_name ]
        request_results_list = []
        for gene_name in gene_name_list:
            try:
                url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28gene:'
                url += gene_name.upper() + '%29%20AND%20%28reviewed%3Atrue%29'
                all_fastas = requests.get(url).text
                fasta_list = re.split(r'\n(?=>)', all_fastas)
                for fasta in fasta_list:
                    fasta_split = re.split('\|', fasta)
                    if 'Homo sapiens' in fasta_split[2]:
                        request_results_list.append(fasta_split[1])
            except Exception:
                print('Uniprot request failed to find gene name ' + gene_name)
        return request_results_list
    
    #
    def map_uniprot_accession_to_gene_name(self, uniprot_accession_list):
        uniprot_accession_chunks = \
            divide_list_into_chunks(uniprot_accession_list, self.chunk_length)
        request_results_list = []
        no_result_list = []
        for x in uniprot_accession_chunks:
            y = self._map_uniprot_accession_to_gene_name(x)
            if y is not None:
                request_results_list.extend(y)
            else:
                for i in range(len(x)):
                    y = self._map_uniprot_accession_to_gene_name([x[i]])
                    if y is not None:
                        request_results_list.extend(y)
                    else:
                        no_result_list.append(x[i])
        return request_results_list, no_result_list