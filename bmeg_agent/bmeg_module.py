import sys
import os

import logging
from bioagents import Bioagent
from indra.statements import Agent
from .bmeg_agent import BMEGAgent
from indra.sources.trips.processor import TripsProcessor
from kqml import KQMLModule, KQMLPerformative, KQMLList, KQMLString, KQMLToken


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('BMEGA')

class BMEGModule(Bioagent):
    name = 'BMEGA'
    tasks = ['FIND-GENE-MUTATION-DATASET','FIND-MUTATION-FREQUENCY', 'FIND-COMMON-PHENOTYPES-FOR-GENES', 'FIND-DRUGS-FOR-MUTATION-DATASET', 'FIND-VARIANTS-FOR-GENES', 'SHOW-MUTATION-DATA']

    def __init__(self, **kwargs):
        self.BA = BMEGAgent()
        # Call the constructor of KQMLModule
        super(BMEGModule, self).__init__(**kwargs)


    def send_display_oncoprint(self, oncoprint_data):
        content = KQMLList('display-oncoprint')

        content.sets('data', str(oncoprint_data))

        self.tell(content)

    def respond_show_mutation_data(self, content):
        """Response content to show-mutation-data request"""
        gene_arg = content.get('GENE')


        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')
        gene_name = gene_names[0]

        disease_arg = content.get('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_kqml_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = _sanitize_disase_name(disease_names[0])
        disease_abbr = self.BA.get_tcga_abbr(disease_name)
        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        oncoprint_data = self.BA.find_variants_for_genes_cbio(gene_names, disease_abbr, "tcga")

        self.send_display_oncoprint(oncoprint_data)

        reply = KQMLList('SUCCESS')

        reply.sets('oncoprint', 'SUCCESS' if len(oncoprint_data) > 0 else 'FAILURE')



        return reply

    def respond_find_mutation_frequency(self, content):
        """Response content to find-mutation-frequency request"""
        gene_arg = content.get('GENE')

        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(gene_arg)

        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')
        gene_name = gene_names[0]

        disease_arg = content.get('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_kqml_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = _sanitize_disase_name(disease_names[0])
        disease_abbr = self.BA.get_tcga_abbr(disease_name)
        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        result = self.BA.find_mutation_frequency(gene_name, disease_abbr)

        if not result:
            return self.make_failure('MISSING_MECHANISM')

        reply = KQMLList('SUCCESS')
        reply.sets('mutfreq', result)

        oncoprint_data = self.BA.find_variants_for_genes_cbio(gene_names, disease_abbr, "tcga")

        self.send_display_oncoprint(oncoprint_data)


        return reply

    def respond_find_variants_for_genes(self, content):
        """Response content to find-variants-for-genes"""
        gene_arg = content.get('GENES')

        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        disease_arg = content.get('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_kqml_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = _sanitize_disase_name(disease_names[0])
        disease_abbr = self.BA.get_tcga_abbr(disease_name)

        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        result = self.BA.find_variants_for_genes_cbio(gene_names, disease_abbr, "tcga")

        if not result:
            return self.make_failure('MISSING_MECHANISM')

        reply = KQMLList('SUCCESS')
        reply.sets('variants', str(result))

        return reply



    def respond_find_drugs_for_mutation_dataset(self, content):
        genes_arg = content.get('GENES')

        if not genes_arg:
            return self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(genes_arg)

        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        dataset_arg = content.gets('DATASET')

        if not dataset_arg:
            dataset_arg = "CCLE"  # default cell line

        result = self.BA.find_drugs_for_mutation_dataset(gene_names, dataset_arg)

        if not result:
            return self.make_failure('NO_DRUGS_FOUND')

        reply = KQMLList('SUCCESS')

        drugs = KQMLList()
        for r in result:
            drugs.append(r)

        reply.set('drugs', drugs)

        return reply


def _sanitize_disase_name(name):
    """Given a disease name returns the sanitized version of it"""
    sanitized_name = name.replace("-", " ").lower()
    return sanitized_name

def _get_kqml_names(kqmlList):
    """Given a kqml list returns the names of sublists in the list"""
    if not kqmlList:
        return None

    arr = kqmlList.data;
    if len(arr) == 0:
        return []

    if not isinstance(arr[0], KQMLList):
        arr = [kqmlList]

    res = list(map(lambda kl: kl.get('NAME').string_value(), arr))

    return res


if __name__ == "__main__":
    BMEGModule(argv=sys.argv[1:])
