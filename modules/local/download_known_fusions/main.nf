process DOWNLOAD_KNOWN_FUSIONS {

   output:
        file("ChimerKB4.xlsx")

   script:
   """
curl -o ChimerKB4.xlsx https://www.kobic.re.kr/chimerdb/downloads?name=ChimerKB4.xlsx
   """
}
