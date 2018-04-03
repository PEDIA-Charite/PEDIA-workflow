from pprint import pprint
import requests
from lib.api import omim
from lib.model import config


def main():
    '''Query list of phenotypic omim ids'''

    config_data = config.ConfigManager("config.ini")
    omim_internal = omim.Omim(config=config_data)

    r = omim_internal.construct_phenotypic_series_mapping()
    print(r)


if __name__ == "__main__":
    main()
