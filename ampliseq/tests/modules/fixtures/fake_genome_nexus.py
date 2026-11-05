#!/usr/bin/env python3
"""
Minimal stand-in for Genome Nexus POST /annotation/genomic, for offline ANNOTATE_GNOMAD tests.
Usage: fake_genome_nexus.py <port> [--lifetime 300]

  17,7579472,7579472,G,C   first request: no gnomAD block (a dropped answer); later: AF + AC/AN,
                           exome nfe 0.74 and genome afr 0.8; ClinVar Benign
  12,25368462,25368462,C,T first request: AC/AN only (exome eas 18394/18394), no ClinVar block;
                           later: ClinVar Pathogenic but no gnomAD block (two half answers)
  anything else            dbSNP-free annotation without gnomAD or ClinVar
The server exits on its own after --lifetime seconds.
"""
import argparse
import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

SEEN = {}


def gnomad(ac_e, an_e, ac_g, an_g, with_af=True, exome_pops=None, genome_pops=None):
    """exome_pops / genome_pops: {'nfe': af} with with_af, else {'eas': (ac, an)}."""
    def source(ac, an, pops):
        data = {'alleleCount': {'ac': ac}, 'alleleNumber': {'an': an}}
        if with_af:
            data['alleleFrequency'] = {'af': round(ac / an, 6)}
        for pop, value in (pops or {}).items():
            if with_af:
                data['alleleFrequency'][f'af_{pop}'] = value
            else:
                data['alleleCount'][f'ac_{pop}'], data['alleleNumber'][f'an_{pop}'] = value
        return data
    return {'gnomadExome': source(ac_e, an_e, exome_pops), 'gnomadGenome': source(ac_g, an_g, genome_pops)}


def clinvar(significance, conflicting=''):
    return {'annotation': {'clinicalSignificance': significance, 'conflictingClinicalSignificance': conflicting}}


def annotate(key):
    SEEN[key] = SEEN.get(key, 0) + 1
    info = {'variant': key}
    answer = {'originalVariantQuery': key, 'successfully_annotated': True}
    if key == '17,7579472,7579472,G,C':
        info['dbsnp'] = {'rsid': 'rs1042522'}
        answer['clinvar'] = clinvar('Benign')
        if SEEN[key] > 1:
            info.update(gnomad(167398, 250534, 19434, 31312, exome_pops={'nfe': 0.74}, genome_pops={'afr': 0.8}))
    elif key == '12,25368462,25368462,C,T':
        if SEEN[key] == 1:
            info.update(gnomad(251165, 251262, 31347, 31394, with_af=False, exome_pops={'eas': (18394, 18394)}))
        else:
            answer['clinvar'] = clinvar('Pathogenic')
    answer['my_variant_info'] = {'annotation': info}
    return answer


class Handler(BaseHTTPRequestHandler):
    def do_POST(self):
        body = json.loads(self.rfile.read(int(self.headers['Content-Length'])))
        keys = [','.join(str(v[f]) for f in ('chromosome', 'start', 'end', 'referenceAllele', 'variantAllele'))
                for v in body]
        payload = json.dumps([annotate(k) for k in keys]).encode()
        self.send_response(200)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def log_message(self, *args):
        pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('port', type=int)
    parser.add_argument('--lifetime', type=float, default=300)
    args = parser.parse_args()
    server = ThreadingHTTPServer(('127.0.0.1', args.port), Handler)
    threading.Timer(args.lifetime, server.shutdown).start()
    server.serve_forever()


if __name__ == '__main__':
    main()
