#! /usr/bin/env python3

import vcf
import httplib2
import json

__author__ = 'ARNOLD Vivienne'

## Aim of this assignment is to annotate the variants with various attributes
## We will use the API provided by "myvariant.info" - more information here: https://docs.myvariant.info
## NOTE NOTE! - check here for hg38 - https://myvariant.info/faq
## 1) Annotate the first 900 variants in the VCF file
## 2) Store the result in a data structure (not in a database)
## 3) Use the data structure to answer the questions
## 4) View the VCF in a browser

class Assignment3:
    
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s\n" % vcf.VERSION)

        ## Call annotate_vcf_file here
        self.vcf_path = "chr16.vcf"

        print("Start investigation of <"+self.vcf_path+">")

        self.annotation_result = self.annotate_vcf_file()

        self.list_of_genes = self.get_list_of_genes()

        self.number_of_variants_modifier = self.get_number_of_variants_modifier()

        self.number_of_variants_with_mutationtaster_annotation = self.get_number_of_variants_with_muttaster_annotation()

        self.number_of_variants_non_synonymous = self.get_number_of_variants_non_synonymous()

        self.url = self.view_vcf_in_browser()

        print("Done")

    def annotate_vcf_file(self):
        '''
        - Annotate the VCF file and store the result in a data structure
        :return: annotation_result_json
        '''    

        ## Build the connection
        h = httplib2.Http()
        headers = {'content-type': 'application/x-www-form-urlencoded'}

        params_pos = []  # List of variant positions
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for counter, record in enumerate(vcf_reader):
                params_pos.append(record.CHROM + ":g." + str(record.POS) + record.REF + ">" + str(record.ALT[0]))
                
                if counter >= 899:
                    break
        
        ## Build the parameters using the list we just built
        params = 'ids=' + ",".join(params_pos) + '&hg38=true'
        
        ## Perform annotation
        res, con = h.request('http://myvariant.info/v1/variant', 'POST', params, headers=headers)
        annotation_result = con.decode('utf-8')

        annotation_result_json = json.loads(annotation_result)

        return (annotation_result_json)

    def get_list_of_genes(self):
        '''
        Print the name of genes in the annotation data set
        :return:
        '''

        count = 0
        gene_list = set()

        for i in self.annotation_result:
            if len(i) > 2:
                if 'cadd' in i:
                    if type(i['cadd']['gene']) is list:
                        count += 1
                        gene_list.add(i['cadd']['gene'][1]['genename'])
                    elif type(i['cadd']['gene']) is dict:
                        count += 1
                        gene_list.add(i['cadd']['gene']['genename'])
                    else:
                        pass

                elif 'dbsnp' in i:
                    if type(i['snpeff']['ann']) is dict:
                        gene_list.add(i['snpeff']['ann']['genename'])
                        count += 1
                    elif type(i['snpeff']['ann']) is list:
                        gene_list.add(i['snpeff']['ann'][0]['genename'])
                        count += 1
                    else:
                        pass

        return sorted(gene_list)

    def get_number_of_variants_modifier(self):
        '''
        Print the number of variants with putative_impact "MODIFIER"
        '''
        number_of_variants_modifier = 0

        for i in self.annotation_result:
            if len(i) > 2:
                if str(i).find('putative_impact\': \'MODIFIER'):
                    number_of_variants_modifier += 1

        return number_of_variants_modifier


    def get_number_of_variants_with_muttaster_annotation(self):
        '''
        Print the number of variants with a 'mutationtaster' annotation
        '''
        variants_with_mutationtaster_annotation = 0

        for i in self.annotation_result:
            if 'dbnsfp' in i:
                if 'mutationtaster' in i['dbnsfp']:
                    variants_with_mutationtaster_annotation += 1
        return variants_with_mutationtaster_annotation

    def get_number_of_variants_non_synonymous(self):
        '''
        Print the number of variants with 'consequence' 'NON_SYNONYMOUS'
        '''
        variants_non_synonymous = 0
        for i in self.annotation_result:
            if 'cadd' in i:
                key, value = "consequence", "NON_SYNONYMOUS"
                if key in i['cadd'] and value == i['cadd']['consequence']:
                    variants_non_synonymous += 1

        return variants_non_synonymous

    def view_vcf_in_browser(self):
        '''
        - Open a browser and go to https://vcf.iobio.io/
        - Upload the VCF file and investigate the details
        :return:
        '''
        ## Document the final URL here
        url = "https://vcf.iobio.io/?species=Human&build=GRCh38"

        return url

    def print_summary(self):
        print("\n"+"*"*80+"\n")

        print("List of Genes:")
        for i in self.list_of_genes:
            print(" * "+i)
        print("")

        print("Frequency of:")
        print(" * variants modifier: "+str(self.number_of_variants_modifier))
        print(" * variants with mutationtaster annotation: "+str(self.number_of_variants_with_mutationtaster_annotation))
        print(" * variants non synonymous: "+str(self.number_of_variants_non_synonymous))

        print("")

        print("For investigation use URL: "+self.url)

        print("\n" + "*" * 80 + "\n")


def main():
    print("Assignment 3\n")
    assignment3 = Assignment3()
    assignment3.print_summary()
    print("Done with assignment 3")


if __name__ == '__main__':
    main()

"""DDX11L10
HBA1
HBA1
LUC7L
NPRL3
RHBDF1
RHBDF1
SNRNP25
Z84721.4"""