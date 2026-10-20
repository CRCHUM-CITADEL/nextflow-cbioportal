#!/usr/bin/env python3
"""
Minimal stand-in for Genome Nexus POST /annotation/genomic, for offline ANNOTATE_GNOMAD tests.
Usage: fake_genome_nexus.py <port> [--lifetime 300]

  17,7579472,7579472,G,C   first request: no gnomAD block (a dropped answer); later: AF + AC/AN
  12,25368462,25368462,C,T AC/AN only, no alleleFrequency
  anything else            dbSNP-free annotation without gnomAD
The server exits on its own after --lifetime seconds.
"""
import argparse
import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

SEEN = {}


def gnomad(ac_e, an_e, ac_g, an_g, with_af=True):
    def source(ac, an):
        data = {'alleleCount': {'ac': ac}, 'alleleNumber': {'an': an}}
        if with_af:
            data['alleleFrequency'] = {'af': round(ac / an, 6)}
        return data
    return {'gnomadExome': source(ac_e, an_e), 'gnomadGenome': source(ac_g, an_g)}


def annotate(key):
    SEEN[key] = SEEN.get(key, 0) + 1
    info = {'variant': key}
    if key == '17,7579472,7579472,G,C':
        info['dbsnp'] = {'rsid': 'rs1042522'}
        if SEEN[key] > 1:
            info.update(gnomad(167398, 250534, 19434, 31312))
    elif key == '12,25368462,25368462,C,T':
        info.update(gnomad(251165, 251262, 31347, 31394, with_af=False))
    return {'originalVariantQuery': key, 'successfully_annotated': True,
            'my_variant_info': {'annotation': info}}


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
