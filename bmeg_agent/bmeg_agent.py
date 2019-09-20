import re
import gripql

import itertools
from scipy import stats
import pdb

import os
# for cbio portal
import http.client, urllib.parse

credentials_dir = os.path.dirname(os.path.realpath(__file__)) + '/../credentials/'

disease_names = {
'ADRENOCORTICAL CARCINOMA':	'ACC',
'BLADDER UROTHELIAL CARCINOMA':	'BLCA',
'BREAST INVASIVE CARCINOMA':	'BRCA',
'INVASIVE BREAST CARCINOMA':	'BRCA',
'BREAST CANCER':	'BRCA',
'BREAST CARCINOMA':	'BRCA',
'CERVICAL CANCER':	'CESC',
'CERVICAL SQUAMOUS CELL CARCINOMA AND ENDOCERVICAL ADENOCARCINOMA':	'CESC',
'ENDOCERVICAL CANCER':	'CESC',
'CHOLANGIOCARCINOMA':	'CHOL',
'COLON ADENOCARCINOMA':	'COAD',
'COLORECTAL CANCER':	'COADREAD',
'LYMPHOID NEOPLASM DIFFUSE LARGE B-CELL LYMPHOMA':	'DLBC',
'GLIOBLASTOMA':	'GBM',
'GLIOBLASTOMA MULTIFORME':	'GBM',
'GLIOBLASTOMA MULTIFORME/BRAIN LOWER GRADE GLIOMA':	'GBMLGG',
'HEAD AND NECK SQUAMOUS CELL CARCINOMA':	'HNSC',
'KIDNEY CHROMOPHOBE':	'KICH',
'PAN-KIDNEY':	'KIPAN',
'KIDNEY RENAL CLEAR CELL CARCINOMA':	'KIRC',
'KIDNEY RENAL PAPILLARY CELL CARCINOMA':	'KIRP',
'KIDNEY CANCER':	'KIRP',
'ACUTE MYELOID LEUKEMIA':	'LAML',
'LEUKEMIA':	'LAML',
'BRAIN LOWER GRADE GLIOMA':'	LGG',
'BRAIN CANCER':	'LGG',
'LIVER HEPATOCELLULAR CARCINOMA':	'LIHC',
'LIVER CANCER':	'LIHC',
'LUNG ADENOCARCINOMA':	'LUAD',
'LUNG CANCER':	'LUAD',
'LUNG SQUAMOUS CELL CARCINOMA':	'LUSC',
'OVARIAN SEROUS CYSTADENOCARCINOMA':	'OV',
'OVARIAN CANCER':	'OV',
'OVARY CANCER':	'OV',
'OVARY SEROUS CYSTADENOCARCINOMA':	'OV',
'PANCREATIC ADENOCARCINOMA':	'PAAD',
'PANCREATIC CANCER':	'PAAD',
'PANCREAS CANCER':	'PAAD',
'PANCREAS ADENOCARCINOMA':	'PAAD',
'PHEOCHROMOCYTOMA AND PARAGANGLIOMA':	'PCPG',
'PROSTATE ADENOCARCINOMA':	'PRAD',
'PROSTATE CANCER':	'PRAD',
'PROSTATIC ADENOCARCINOMA':	'PRAD',
'PROSTATIC CANCER':	'PRAD',
'RECTUM ADENOCARCINOMA':	'READ',
'RECTAL ADENOCARCINOMA':	'READ',
'RECTAL CANCER':	'READ',
'RECTUM CANCER':	'READ',
'SARCOMA':	'SARC',
'SKIN CUTANEOUS MELANOMA':	'SKCM',
'MELANOMA':	'SKCM',
'STOMACH ADENOCARCINOMA':	'STAD',
'STOMACH CANCER':	'STAD',
'ESOPHAGUS CARCINOMA':	'STES',
'ESOPHAGUS CANCER':	'STES',
'ESOPHAGEAL CANCER':	'STES',
'STOMACH AND ESOPHAGEAL CARCINOMA':	'STES',
'ESOPHAGEAL CARCINOMA':	'STES',
'TESTIS CANCER':	'TGCT',
'TESTICULAR GERM CELL TUMORS':	'TGCT',
'TESTICULAR CANCER':	'TGCT',
'THYROID CANCER':	'THCA',
'THYROID CARCINOMA':	'THCA',
'UTERINE CORPUS ENDOMETRIAL CARCINOMA':	'UCEC',
'UTERINE CORPUS ENDOMETRIAL CANCER':	'UCEC',
'UTERINE CARCINOSARCOMA':	'UCS',
'UTERINE CANCER':	'UCS',
'UVEAL MELANOMA':	'UVM',
}





class BMEGAgent:
    def __init__(self):
        conn = gripql.Connection('https://bmeg.io/api', credential_file= credentials_dir+ "bmeg_credentials.json")
        self.O = conn.graph("bmeg_rc2")
        print("Connected to bmeg")

    def get_tcga_abbr(self, long_name):
        return disease_names[long_name.upper()]

    def find_mutation_frequency(self, gene, disease):

        # find gene id


        for i in self.O.query().V().hasLabel("Gene").has(gripql.eq("symbol", gene)):
            gene_id = i.gid

        # Find all the tumor aliquots in tcga for the disease
        proj_id = "Project:TCGA-" + disease.upper()
        q = self.O.query().V(proj_id).out("cases").out("samples").out("aliquots").has(
            gripql.eq("gdc_attributes.sample_type", "Primary Tumor")).as_("sample").out("somatic_callsets").select("sample")
        all_aliquots = []
        for row in q:
            all_aliquots.append(row.gid)
        if len(all_aliquots) == 0:
            return 0



        q = self.O.query().V(all_aliquots).as_("sample").out("somatic_callsets").outE("alleles")
        q = q.has(gripql.eq("ensembl_gene", gene_id)).as_("variant")
        q = q.distinct("_from")

        mut_samples = []
        for i in q:
            mut_samples.append(i.gid)

        freq = (float(len(mut_samples)) / float(len(all_aliquots))) #* 100

        # print (freq)
        return freq


    def find_drugs_for_mutation_dataset(self, genes, dataset):
        dataset = dataset.upper()
        q = self.O.query().V().has(getGripqlProjIdRange(dataset)).out("cases").distinct("_gid")

        all_cases = []
        for row in q:
            all_cases.append(row.gid)

        gene_ids = {}
        for i in self.O.query().V().hasLabel("Gene").has(gripql.within("symbol", genes)):
            gene_ids[i.data.symbol] = i.gid

        mut_cases = {}
        norm_cases = {}

        q = self.O.query().V(all_cases).as_("ds")

        if dataset != "CCLE":
            q = q.out("same_as").has(getGripqlProjIdRange("Project:CCLE"))

        q = q.out("samples").out("aliquots").out("somatic_callsets")
        q = q.outE("alleles").has(gripql.within("ensembl_gene", list(gene_ids.values())))
        q = q.render({"case" : "$ds._gid", "gene" : "$._data.ensembl_gene"})

        for res in q:
            mut_cases[res.gene] = mut_cases.get(res.gene, set()) | set([res.case])

        #get CCLE samples without mutation
        for i in gene_ids.values():
            norm_cases[i] = list(set(all_cases).difference(mut_cases[i]))

            print( "%s Positive Set: %d" % (i, len(mut_cases[i])) )
            print( "%s Negative Set: %d" % (i, len(norm_cases[i])) )

        names = {}

        area_metric = "act_area" if dataset == "CCLE" else "auc"

        pos_response = {}
        for g in gene_ids.values():
            pos_response[g] = {}
            q = self.O.query().V(list(mut_cases[g])).as_("a").out("samples").out("aliquots")
            q = q.out("drug_response").as_("a").out("compounds").as_("b")
            q = q.select(["a", "b"])
            for row in q:
                if hasattr(row["a"]["data"], area_metric):
                    v = row["a"]["data"][area_metric]
                else:
                    v = 0

                id = row["b"]["gid"]
                names[id] = row["b"]["data"]["name"]
                if id not in pos_response[g]:
                    pos_response[g][id] = [ v ]
                else:
                    pos_response[g][id].append(v)

        neg_response = {}
        for g in gene_ids.values():
            neg_response[g] = {}
            q = self.O.query().V(list(norm_cases[g])).as_("a").out("samples").out("aliquots")
            q = q.out("drug_response").as_("a").out("compounds").as_("b")
            q = q.select(["a", "b"])
            for row in q:
                if hasattr(row["a"]["data"], area_metric):
                    v = row["a"]["data"][area_metric]
                else:
                    v = 0

                id = row["b"]["gid"]
                names[id] = row["b"]["data"]["name"]
                if id not in neg_response[g]:
                    neg_response[g][id] = [ v ]
                else:
                    neg_response[g][id].append(v)


        drugs = set(itertools.chain.from_iterable( i.keys() for i in pos_response.values() ))
        out = []
        for drug in drugs:
            for g in gene_ids.values():
                if drug in pos_response[g] and drug in neg_response[g]:
                    row = {"drug" : drug, "mutation" : g}
                    mut_values = pos_response[g][drug]
                    norm_values = neg_response[g][drug]
                    if len(mut_values) > 5 and len(norm_values) > 5:
                        s = stats.ttest_ind(mut_values, norm_values, equal_var=False)
                        if s.pvalue <= 0.05 and s.statistic > 0:  # means drug is significantly effective
                            out.append(names[drug])
        return out

    def find_variants_for_genes_cbio(self, genes, disease, dataset= "tcga"):
        """

        :param genes:
        :param disease: TCGA abbreviation of the study
        :param dataset:
        :return:
        """

        disease = disease.lower()
        headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/plain"}

        gene_list = ','.join(genes)

        genetic_profile_id = disease + "_" + dataset + "_mutations"
        params = urllib.parse.urlencode({'cmd': 'getMutationData', 'gene_list': gene_list, 'genetic_profile_id': genetic_profile_id})


        conn = http.client.HTTPConnection("www.cbioportal.org")

        conn.request("POST", "/webservice.do?", params, headers)



        r2 = conn.getresponse()



        if r2.status != 200:
            return -1

        data = r2.read().decode("utf-8")



        # Now get CNA data

        case_set_id = disease + "_tcga_all"
        genetic_profile_id =  disease + "_" + dataset + "_gistic"
        params = urllib.parse.urlencode(
            {'cmd': 'getProfileData', 'gene_list': gene_list, 'genetic_profile_id': genetic_profile_id,
             'case_set_id': case_set_id})

        conn.request("POST", "/webservice.do?", params, headers)

        r2 = conn.getresponse()

        if r2.status != 200:
            return format_cbio_mutation_output(data)


        data2 = r2.read().decode("utf-8")

        conn.close()

        out_dict = dict()

        format_cbio_mutation_output(data, out_dict)

        # uses old out_dict to append new data
        out = format_cbio_cna_output(data2, out_dict)



        return out

    # Looks like this method is never used so just commented it out for now.
    # The changes that are done for adopting 'bmeg_rc2' is not applied to this method.

    # def find_variants_for_genes_bmeg(self, genes, disease, dataset= "tcga"):
    #
    #     gene_ids = {}
    #     for g in genes:
    #         for i in self.O.query().V().hasLabel("Gene").has(gripql.eq("symbol", g)):
    #             gene_ids[g] = i.gid
    #
    #
    #
    #         # Find all the tumor aliquots in tcga for the disease
    #     proj_id = "Project:TCGA-" + disease
    #     q = self.O.query().V(proj_id).in_("InProject").in_("SampleFor").in_("AliquotFor").has(
    #         gripql.eq("gdc_attributes.sample_type", "Primary Tumor")).as_("sample").in_("CallsetFor").select("sample")
    #     all_aliquots = []
    #     for row in q:
    #         all_aliquots.append(row.gid)
    #     if len(all_aliquots) == 0:
    #         return 0
    #
    #     q = self.O.query().V(all_aliquots).as_("sample").in_("CallsetFor").outV("AlleleCall")
    #
    #     #
    #     # q = q.has(gripql.within("ensembl_gene", list(gene_ids.values()))).as_("variant")
    #     #
    #     # q = q.distinct("gid")
    #
    #     mutations = []
    #     for i in q:
    #         if hasattr(i['data'], 'effect'): # not all rows have 'effect' attribute
    #             mutations.append(i.data.effect)
    #
    #
    #     return mutations


def format_cbio_cna_output(data, out_dict):
    """

        :param data: Is tab separated
        :return:
        """


    lines = data.split("\n")

    samples = lines[2].split("\t")[2:]

    for line in lines[3:]:

        words = line.split("\t")

        if len(words) < 2:
            continue

        gene = words[1]

        for i, s in enumerate(samples):
            data = {}
            data['sample'] = s

            data['disp_cna'] = map_to_oncoprint_cna(words[i+2])


            if gene in out_dict:
                out_dict[gene]['data'].append(data)
            else:
                out_dict[gene] = {}
                out_dict[gene]['data'] = [data]
                out_dict[gene]['gene'] = gene
                out_dict[gene]['desc'] = "Annotation for " + gene


    out = []
    for gene in out_dict:
        out.append(out_dict[gene])



    return out


def format_cbio_mutation_output( data ,out_dict):
    """

    :param data: Is tab separated
    :return:
    """
    # print(data)



    lines = data.split("\n")

    for line in lines[2:]:
        words = line.split("\t")
        if len(words) < 6:
            continue
        gene = words[1]
        data = {'sample':'', 'disp_mut':''}
        data['sample'] = words[2]

        data['disp_mut'] = map_to_oncoprint_mutation(words[5])

        if gene in out_dict:
            out_dict[gene]['data'].append(data)
        else:
            out_dict[gene] = {}
            out_dict[gene]['data'] = [data]
            out_dict[gene]['gene'] = gene
            out_dict[gene]['desc'] = "Annotation for " + gene

    out = []
    for gene in out_dict:
        out.append(out_dict[gene])


    return out


def map_to_oncoprint_cna(ind_str):
    if ind_str.lower() == 'nan':
        return '*'
    ind = int(ind_str)

    cna_type = '*'
    if ind == -2:
        cna_type = "homdel"
    elif ind == -1:
        cna_type = "hetloss"
    elif ind == 0:
        cna_type = "diploid"
    elif ind == 1:
        cna_type = "gain"
    elif ind == 2:
        cna_type = "amp"

    return cna_type


def map_to_oncoprint_mutation(cbio_mut):
    cbio_mut = cbio_mut.lower()

    oncoprint_mut = '*'
    if 'missense' in cbio_mut:
        oncoprint_mut = 'missense'
    elif 'inframe' in cbio_mut:
        oncoprint_mut = 'inframe'
    elif 'promote' in cbio_mut:
        oncoprint_mut = 'promoter'
    elif 'frame_shift' in cbio_mut:
        oncoprint_mut = 'trunc'
    elif 'nonsense' in cbio_mut:
        oncoprint_mut = 'trunc'
    elif 'splice' in cbio_mut:
        oncoprint_mut = 'trunc'
    elif 'amp' in cbio_mut:
        oncoprint_mut = 'amp'

    return oncoprint_mut

def getGripqlProjIdRange(prefix):
    """
    Returns a gripql range where project_id starts with given prefix assuming that
    the following part of project_id will start by "_" chacter that is followed by
    a letter.
    """
    return gripql.between("project_id", prefix + "_A", prefix + "_zz")

# ba = BMEGAgent()
# ba.find_mutation_frequency("TP53", "OV")

# print(ba.find_variants_for_genes_cbio(['PTEN'],'BRCA','tcga'))
# print(ba.find_variants_for_genes_cbio(['BRAF'],'OV','tcga'))
# ba.find_variants_for_genes_cbio(['EGFR', 'PTEN'],'glioblastoma','tcga')
# ba.find_variants_for_genes_cbio(['BRCA1, CCNE1'],'OV','tcga')

# ba.find_mutations_on_gene("BRAF")
# ba.find_drugs_for_mutation_dataset(["CDKN2A", "PTEN", "TP53",], "CTRP")
# ba.find_drugs_for_mutation_dataset(["TP53"], "CCLE")


# for bmeg

# ba.find_variants_for_genes_bmeg(['TP53'],'OV','tcga')
