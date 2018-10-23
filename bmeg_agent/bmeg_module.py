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
    tasks = ['FIND-GENE-MUTATION-DATASET', ]

    def __init__(self, **kwargs):
        self.BA = BMEGAgent()
        # Call the constructor of KQMLModule
        super(BMEGModule, self).__init__(**kwargs)



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
