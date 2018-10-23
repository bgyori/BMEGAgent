import json
from kqml import KQMLList, KQMLString, KQMLPerformative
from indra.statements import stmts_from_json

from bmeg_agent.bmeg_agent import BMEGAgent
from bmeg_agent.bmeg_module import BMEGModule
from bioagents.tests.integration import _IntegrationTest
from bioagents.tests.util import ekb_kstring_from_text, ekb_from_text, get_request
import time

ca = BMEGAgent.BMEGA()


class TestBMEGAgent(_IntegrationTest):
    def __init__(self, *args):
        super(TestBMEGAgent, self).__init__(BMEGModule)


class TestMutFreq(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutFreq, self).__init__(BMEGModule)

    def create_message_OV(self):
        content = KQMLList('FIND-MUTATION-FREQUENCY')
        gene = ekb_kstring_from_text('TP53')
        disease = ekb_from_text('Ovarian serous cystadenocarcinoma')
        content.set('gene', gene)
        content.set('disease', disease)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_OV(self, output):
        assert output.head() == 'SUCCESS', output
        mut_freq = output.gets('mutfreq')
        assert mut_freq.startswith('0.81')
