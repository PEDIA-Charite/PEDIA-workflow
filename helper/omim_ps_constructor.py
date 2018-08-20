# snipplets from omim api for phenotypic series creation
    def api_query_phenotype_mapping(self, mim_pheno: str):
        '''Map phenotype to phenotypic series number.'''
        entry_list = self.api_query_entry(mim_pheno)
        phenotypic_numbers = [
            v['phenotypeMap']['phenotypicSeriesNumber']
            for e in entry_list
            if 'phenotypeMapList' in e
            for v in e['phenotypeMapList']
            if 'phenotypicSeriesExists' in e and e['phenotypicSeriesExists']
        ]
        return phenotypic_numbers

    def construct_phenotypic_series_mapping(self):
        '''Construct a dict with mapping of phenotypic omim ids to
        phenotypic series ids.'''

        output_path = os.path.join(self._mimdir, "mim_to_ps.json")
        if os.path.exists(output_path):
            with open(output_path, "r") as old_json:
                old = json.load(old_json)
                if old['morbidmap_hash'] == self.morbidmap_hash:
                    LOGGER.warning(
                        ("Phenotypic series mapping with same base "
                         "morbidmap has already been constructed.")
                    )
                    return
        mim_numbers = self.morbidmap[
            "phen_mim_number"].dropna().drop_duplicates()

        mim_ps_mapping = {
            num: self.api_query_phenotype_mapping(num)
            for num in mim_numbers
        }
        mim_ps_data = {
            'morbidmap_hash': self.morbidmap_hash,
            'mapping': mim_ps_mapping
        }
        with open(output_path, "w") as new_json:
            json.dump(mim_ps_data, new_json)
        return mim_ps_data
