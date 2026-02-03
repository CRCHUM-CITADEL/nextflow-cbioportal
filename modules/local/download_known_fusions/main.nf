process DOWNLOAD_KNOWN_FUSIONS {
   
   output:
        file("ChimerKB4.xlsx")

   script:
   """
curl https://www.kobic.re.kr/chimerdb/downloads?name=ChimerKB4.xlsx --name ChimerKB4.xlsx
   """
}
