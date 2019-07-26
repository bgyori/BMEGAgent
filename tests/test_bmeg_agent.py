import json
from kqml import KQMLList, KQMLString, KQMLPerformative
from indra.statements import stmts_from_json

from bmeg_agent.bmeg_agent import BMEGAgent
from bmeg_agent.bmeg_module import BMEGModule
from bioagents.tests.integration import _IntegrationTest
from bioagents.tests.util import ekb_kstring_from_text, ekb_from_text, get_request, agent_clj_from_text

class TestShowMutData(_IntegrationTest):
    def __init__(self, *args):
        super(TestShowMutData, self).__init__(BMEGModule)

    def create_message_1(self):
        content = KQMLList('SHOW-MUTATION-DATA')
        gene = agent_clj_from_text('TP53')
        disease = agent_clj_from_text('Ovarian serous cystadenocarcinoma')
        content.set('gene', gene)
        content.set('disease', disease)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_1(self, output):
        assert output.head() == 'SUCCESS', output
        oncoprint = output.gets('oncoprint')
        assert oncoprint == 'SUCCESS'


class TestMutFreq(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutFreq, self).__init__(BMEGModule)

    def create_message_OV(self):
        content = KQMLList('FIND-MUTATION-FREQUENCY')
        gene = agent_clj_from_text('TP53')
        disease = agent_clj_from_text('Ovarian serous cystadenocarcinoma')
        content.set('gene', gene)
        content.set('disease', disease)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_OV(self, output):
        assert output.head() == 'SUCCESS', output
        mut_freq = output.gets('mutfreq')
        assert mut_freq.startswith('0.9')


class TestDrugMutationDataset(_IntegrationTest):
    def __init__(self, *args):
        super(TestDrugMutationDataset, self).__init__(BMEGModule)
        self.timeout = 200

    def create_message_1(self):
        content = KQMLList('FIND-DRUGS-FOR-MUTATION-DATASET')
        genes = agent_clj_from_text('TP53')
        content.set('genes', genes)
        content.set('dataset', "CTRP")

        msg = get_request(content)
        return msg, content

    def check_response_to_message_1(self, output):
        assert output.head() == 'SUCCESS', output

        drugs = output.get('drugs')

        assert  'selumetinib' in drugs

class TestVariantsCbioportal(_IntegrationTest):
    def __init__(self, *args):
        super(TestVariantsCbioportal, self).__init__(BMEGModule)

    def create_message_1(self):
        content = KQMLList('FIND-VARIANTS-FOR-GENES')
        genes = agent_clj_from_text('EGFR and PTEN')
        content.set('genes', genes)

        disease = agent_clj_from_text('glioblastoma')
        content.set('disease', disease)
        content.sets('dataset', "tcga")

        msg = get_request(content)
        return msg, content

    def check_response_to_message_1(self, output):
        assert output.head() == 'SUCCESS', output

        variants = output.gets('variants')

        assert 'missense' in variants

    # def create_message_2(self):
    #     content = KQMLList('FIND-VARIANTS-FOR-GENES')
    #     genes = ekb_from_text('CCNE1')
    #     content.sets('genes', str(genes))
    #
    #     disease = ekb_from_text('ovarian cancer')
    #     content.sets('disease', disease)
    #     content.sets('dataset', "tcga")
    #
    #
    #     msg = get_request(content)
    #     return msg, content
    #
    # def check_response_to_message_2(self, output):
    #     assert output.head() == 'SUCCESS', output
    #
    #     variants = output.gets('variants')
    #
    #     assert 'amplification' in variants


