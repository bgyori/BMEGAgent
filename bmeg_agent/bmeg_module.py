import sys
import os

import logging
from bioagents import Bioagent
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
        gene_arg = content.gets('GENE')


        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_term_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')
        gene_name = gene_names[0]

        disease_arg = content.gets('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_term_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = disease_names[0].replace("-", " ").lower()
        disease_abbr = self.BA.get_tcga_abbr(disease_name)
        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')


        gene_list = []
        for gene_name in gene_names:
            gene_list.append(str(gene_name))

        oncoprint_data = self.BA.find_variants_for_genes_cbio(gene_list, disease_abbr, "tcga")

        self.send_display_oncoprint(oncoprint_data)

        reply = KQMLList('SUCCESS')

        reply.sets('oncoprint', 'SUCCESS' if len(oncoprint_data) > 0 else 'FAILURE')



        return reply

    def respond_find_mutation_frequency(self, content):
        """Response content to find-mutation-frequency request"""
        gene_arg = content.gets('GENE')

        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_term_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')
        gene_name = gene_names[0]

        disease_arg = content.gets('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_term_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = disease_names[0].replace("-", " ").lower()
        disease_abbr = self.BA.get_tcga_abbr(disease_name)
        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        result = self.BA.find_mutation_frequency(gene_name, disease_abbr)

        if not result:
            return self.make_failure('MISSING_MECHANISM')

        reply = KQMLList('SUCCESS')
        reply.sets('mutfreq', result)

        gene_list = []
        for gene_name in gene_names:
            gene_list.append(str(gene_name))

        oncoprint_data = self.BA.find_variants_for_genes_cbio(gene_list, disease_abbr, "tcga")

        self.send_display_oncoprint(oncoprint_data)


        return reply

    def respond_find_variants_for_genes(self, content):
        """Response content to find-variants-for-genes"""
        gene_arg = content.gets('GENES')

        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_term_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        gene_list = []
        for gene_name in gene_names:
            gene_list.append(str(gene_name))

        disease_arg = content.gets('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_term_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = disease_names[0].replace("-", " ").lower()
        disease_abbr = self.BA.get_tcga_abbr(disease_name)

        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        result = self.BA.find_variants_for_genes_cbio(gene_list, disease_abbr, "tcga")

        if not result:
            return self.make_failure('MISSING_MECHANISM')

        reply = KQMLList('SUCCESS')
        reply.sets('variants', str(result))

        return reply



    def respond_find_drugs_for_mutation_dataset(self, content):
        genes_arg = content.gets('GENES')

        if not genes_arg:
            return self.make_failure('MISSING_MECHANISM')

        gene_names = _get_term_names(genes_arg)

        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        gene_list = []
        for gene_name in gene_names:
            gene_list.append(str(gene_name))


        dataset_arg = content.gets('DATASET')

        if not dataset_arg:
            dataset_arg = "CCLE"  # default cell line

        result = self.BA.find_drugs_for_mutation_dataset(gene_list, dataset_arg.lower())

        if not result:
            return self.make_failure('NO_DRUGS_FOUND')

        reply = KQMLList('SUCCESS')

        drugs = KQMLList()
        for r in result:
            drugs.append(r)

        reply.set('drugs', drugs)

        # drugs = KQMLList.from_string(drugs.to_string())
        # reply.set('drugs', drugs)

        return reply



def _get_term_names(term_str):
    """Given an ekb-xml returns the names of genes in a list"""

    tp = TripsProcessor(term_str)
    terms = tp.tree.findall('TERM')
    if not terms:
        return None

    agent_names = []
    for term in terms:
        term_id = term.attrib['id']
        agent = tp._get_agent_by_id(term_id, None)

        if agent is not None:
            if isinstance(agent, list):
                for a in agent:
                    if a.name:
                        agent_names.append(a.name)
            else:
                agent_names.append(agent.name)

    if len(agent_names) == 0:
        return None

    return agent_names



if __name__ == "__main__":
    BMEGModule(argv=sys.argv[1:])
