import re
import gripql
conn = gripql.Connection('http://bmeg.io')
O = conn.graph("bmeg")
import itertools
from scipy import stats
import pandas


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
        print("Connected to bmeg")

    def get_tcga_abbr(self, long_name):
        return disease_names[long_name.upper()]

    def find_mutation_frequency(self, gene, disease):

        #
        q = O.query().V().where(gripql.eq("_label", "Biosample"))

        q = q.where(gripql.and_(gripql.eq("source", "tcga"), gripql.eq("disease_code", disease))).render({"id": "_gid"})
        all_samples = []
        for row in q:
            all_samples.append(row.id)

        if len(all_samples) == 0:
            return 0
        gene_id = 0
        for i in O.query().V().where(gripql.eq("_label", "Gene")).where(gripql.eq("symbol", gene)):
            gene_id = i.gid

        mut_samples = []

        # get TCGA samples with mutation

        for i in O.query().V(gene_id).in_("variantIn").out("variantCall").out("callSetOf").where(
                gripql.in_("_gid", all_samples)).render({"gid": "_gid"}):
            mut_samples.append(i.gid)
        #

        freq = (float(len(mut_samples)) / float(len(all_samples))) * 100

        return freq
        # return 2


    def find_common_phenotypes_for_genes(self, genes):
        """
        Looks at the mutations on these genes and finds common phenotypes for the mutations
        :param genes: Gene names as a list
        :return:
        """

        gene_ids = {}
        for g in genes:
            for i in O.query().V().where(gripql.eq("_label", "Gene")).where(gripql.eq("symbol", g)):
                gene_ids[g] = i.gid


        phenotypes = []
        ind = 0
        for g, i in gene_ids.items():
            q = O.query().V(i).in_("variantIn").in_("featureOf").out("phenotypeOf")

            phenotypes.append([])

            for k in q:
                if k['data']['description'] not in phenotypes[ind]:
                    phenotypes[ind].append(k['data']['description'])
            ind +=1

        intSet = set(phenotypes[0])
        for i in range(1, len(phenotypes)): #get intersection of two sets
           intSet = intSet & set(phenotypes[i])

        return list(intSet)



    # def find_mutations_on_gene(self, gene):
    #     """
    #     Returns the the mutations on a gene
    #     :param genes: Gene names as a list
    #     :return:
    #     """
    #
    #     gene_id = 0
    #     for i in O.query().V().where(gripql.eq("_label", "Gene")).where(gripql.eq("symbol", gene)):
    #         gene_id = i.gid
    #
    #     #
    #     q = O.query().V(gene_id).in_("variantIn").in_("featureOf").out("environmentFor")
    #     for k in q: #k['data']['description'] in q:
    #         print(k)
    #
    #
    #     #
    #     #     phenotypes.append([])
    #     #
    #     #     for k in q:
    #     #         if k['data']['description'] not in phenotypes[ind]:
    #     #             phenotypes[ind].append(k['data']['description'])
    #     #     ind +=1
    #     #
    #     # intSet = set(phenotypes[0])
    #     # for i in range(1, len(phenotypes)): #get intersection of two sets
    #     #    intSet = intSet & set(phenotypes[i])
    #     #
    #     # return list(intSet)


    def find_drugs_for_mutation_dataset(self, genes, dataset):

        q = O.query().V().where(gripql.eq("_label", "Biosample"))
        q = q.where(gripql.and_(gripql.eq("source", dataset))).render({"id": "_gid"})
        all_samples = []
        for row in q:
            all_samples.append(row.id)

        # GENES = ["CDKN2A", "PTEN", "TP53", "SMAD4"]
        gene_ids = {}
        for g in genes:
            for i in O.query().V().where(gripql.eq("_label", "Gene")).where(gripql.eq("symbol", g)):
                gene_ids[g] = i.gid


        #Scan <dataset> cell lines based on mutation status
        mut_samples = {}
        norm_samples = {}
        for g, i in gene_ids.items():
            # get CCLE samples with mutation
            mut_samples[g] = set(k['gid'] for k in
                                 O.query().V(i).in_("variantIn").out("variantCall").out("callSetOf").where(
                                     gripql.in_("_gid", all_samples)).render({"gid": "_gid"}))

            # get CCLE samples without mutation
            norm_samples[g] = list(set(all_samples).difference(mut_samples[g]))

            print("%s Positive Set: %d" % (g, len(mut_samples[g])))
            print("%s Negative Set: %d" % (g, len(norm_samples[g])))


        # Get response values for the positive set (samples with mutation) and collect AUC value by drug
        pos_response = {}
        compound = {}
        for g in genes:
            pos_response[g] = {}

            for row in O.query().V(list(mut_samples[g])).in_("responseFor").mark("a").out("responseTo").mark("b").select(["a", "b"]):

                for v in row['a']['data']['summary']:
                    if v['type'] == "AMAX":
                        id = row['b']['gid']
                        compound[id] = row['b']['data']['name']

                        if id not in pos_response[g]:
                            pos_response[g][id] = [ v["value"] ]
                        else:
                            pos_response[g][id].append(v["value"])


        #Get response values for the negative set (samples without mutation) and collect AUC value by drug
        neg_response = {}
        for g in genes:
            neg_response[g] = {}
            for row in O.query().V(norm_samples[g]).in_("responseFor").mark("a").out("responseTo").mark("b").select(
                    ["a", "b"]):
                for v in row['a']['data']['summary']:
                    if v['type'] == "AMAX":
                        id = row['b']['gid']
                        compound[id] = row['b']['data']['name']

                        if id not in neg_response[g]:
                            neg_response[g][id] = [v["value"]]
                        else:
                            neg_response[g][id].append(v["value"])

        #Collect t-test statistics
        drugs = set(itertools.chain.from_iterable(i.keys() for i in pos_response.values()))
        out = []
        for drug in drugs:
            for g in genes:
                if drug in pos_response[g] and drug in neg_response[g]:
                    # row = {"drugId": drug, "drugName": compound[drug], "gene": g}

                    mut_values = pos_response[g][drug]
                    norm_values = neg_response[g][drug]
                    if len(mut_values) > 5 and len(norm_values) > 5:
                        s = stats.ttest_ind(mut_values, norm_values, equal_var=False)
                        if s.pvalue <= 0.05 and s.statistic > 0: # means drug is significantly effective
                            out.append(compound[drug])
                            # row["t-statistic"] = s.statistic
                            # row["t-pvalue"] = s.pvalue
                            # out.append(row)



        # print(out)
        return out
        # pd = pandas.DataFrame(out, columns=["drug id", "drug name", "mutation", "t-statistic", "t-pvalue"])


        # print(pd)
        # return pd


# ba = BMEGAgent()
# ba.find_mutations_on_gene("BRAF")
# ba.find_common_phenotypes_for_genes(["BRAF", 'AKT1'])
# ba.find_drugs_for_gene_mutation_dataset(["CDKN2A", "PTEN", "TP53",], "ccle")
# ba.find_drugs_for_gene_mutation_dataset(["CDKN2A", "PTEN", "TP53", "SMAD4"], "ccle")
