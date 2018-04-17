'''Create an omim mapping based on information from Face2Gene using the
Face2Gene library.'''

from pprint import pprint
import json

from lib.api.face2gene import Face2Gene
from lib.model.config import ConfigManager

config_data = ConfigManager()

f2g_session = Face2Gene(config=config_data)

s_list = f2g_session.browse_all_syndromes()

with open("f2g_library_dump.json", "w") as f2g_dump:
    json.dump(s_list, f2g_dump, indent=4)
