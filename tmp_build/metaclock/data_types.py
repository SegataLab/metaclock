class ContigsType(object):
    def __init__(self, reference_genome, samples, **kwargs):
        self.reference_genome = reference_genome
        self.samples = samples

        if 'parameter_set' not in kwargs:
            raise TypeError('parameter_set is not provided')
        param_set = kwargs['parameter_set']

        self.intermediate = kwargs.get('intermediate', '')
        self.parameter_set = {}
        self.param_set['homolog_length'] = param_set.get(
                                            'homolog_length', 500)
        self.param_set['homolog_identity'] = param_set.get(
                                            'homolog_identity', 95.0)
        self.param_set['blastn_threads'] = param_set.get('blastn_threads', 6)


class ReadsType(object):
    def __init__(self, reference_genome, samples):
        self.input_type = 'reads'
        self.reference_genome = reference_genome
        self.samples = samples

        self.validate()

    def validate(self):
        if type(self.reference_genome) is not str:
            raise TypeError("reference_genome should be str type")


class AncientReadsType(ReadsType):
    def __init__(self, reference_genome, samples, **kwargs):
        super().__init__(reference_genome, samples)

        self.age_type = 1
        self.intermediate = kwargs.get('intermediate', 'intermediates')
        self.search_report_mode = kwargs.get('search_report_mode', '-k,5')
        self.bowtie2_threads = kwargs.get('bowtie2_threads', 1)
        self.minimum_mapping_quality = kwargs.get('minimum_mapping_quality', 30)
        self.minimum_mapping_length = kwargs.get('minimum_mapping_length', 30)
        self.maximum_snp_edit_distance = kwargs.get('maximum_snp_edit_distance', 0.03)
        self.nproc = kwargs.get('nproc', 1)
        self.minimum_coverage = kwargs.get('minimum_coverage', 5)
        self.trim_distance = kwargs.get('trim_distance', '5:5')
        self.dominant_allele_frequency = kwargs.get('dominant_allele_frequency', 0.8)
        self.output_trimmed_reads = kwargs.get('output_trimmed_reads', 0)

    def validate(self):
        super().validate()

        if type(self.intermediate) is not str:
            raise TypeError("intermediate should be str type")

        if type(self.search_report_mode) is not str:
            raise TypeError("search_report_mode should be str type")

        if type(self.bowtie2_threads) is not str:
            raise TypeError("bowtie2_threads should be str type")

        if type(self.minimum_mapping_quality) is not int:
            raise TypeError("minimum_mapping_quality should be int type")

        if type(self.minimum_mapping_length) is not int:
            raise TypeError("minimum_mapping_length should be int type")

        if type(self.maximum_snp_edit_distance) is not int:
            raise TypeError("maximum_snp_edit_distance should be float type")

        if type(self.nproc) is not int:
            raise TypeError("nproc should be int type")

        if type(self.minimum_coverage) is not int:
            raise TypeError("minimum_coverage should be int type")

        if type(self.trim_distance) is not str:
            raise TypeError("trim_distance should be str type")

        if not re.match('[0-9]:[0-9]', self.trim_distance):
            raise TypeError("trim_distance should ' \
                            'match re pattern `[0-9]:[0-9]`")

        if type(self.dominant_allele_frequency) is not float \
            and self.dominant_allele_frequency is not 0:
            raise TypeError("dominant_allele_frequency should be float type")

        if type(self.output_trimmed_reads) is not int:
            raise TypeError("output_trimmed_reads should be int type")


class ModernReadsType(ReadsType):
    def __init__(self, reference_genome, samples, **kwargs):
        super().__init__(reference_genome, samples)

        self.age_type = 2
        self.intermediate = kwargs.get('intermediate', 'intermediates')
        self.search_report_mode = kwargs.get('search_report_mode', '-k,5')
        self.bowtie2_threads = kwargs.get('bowtie2_threads', 1)
        self.minimum_mapping_quality = kwargs.get('minimum_mapping_quality', 30)
        self.minimum_mapping_length = kwargs.get('minimum_mapping_length', 30)
        self.maximum_snp_edit_distance = kwargs.get('maximum_snp_edit_distance', 0.03)
        self.nproc = kwargs.get('nproc', 1)
        self.minimum_coverage = kwargs.get('minimum_coverage', 5)
        self.trim_distance = kwargs.get('trim_distance', '5:5')
        self.dominant_allele_frequency = kwargs.get('dominant_allele_frequency', 0.8)

    def validate(self):
        super().validate()

        if type(self.intermediate) is not str:
            raise TypeError("intermediate should be str type")

        if type(self.search_report_mode) is not str:
            raise TypeError("search_report_mode should be str type")

        if type(self.bowtie2_threads) is not str:
            raise TypeError("bowtie2_threads should be str type")

        if type(self.minimum_mapping_quality) is not int:
            raise TypeError("minimum_mapping_quality should be int type")

        if type(self.minimum_mapping_length) is not int:
            raise TypeError("minimum_mapping_length should be int type")

        if type(self.maximum_snp_edit_distance) is not int:
            raise TypeError("maximum_snp_edit_distance should be float type")

        if type(self.nproc) is not int:
            raise TypeError("nproc should be int type")

        if type(self.minimum_coverage) is not int:
            raise TypeError("minimum_coverage should be int type")

        if type(self.trim_distance) is not str:
            raise TypeError("trim_distance should be str type")

        if not re.match('[0-9]:[0-9]', self.trim_distance):
            raise TypeError("trim_distance should ' \
                            'match re pattern `[0-9]:[0-9]`")

        if type(self.dominant_allele_frequency) is not float \
            and self.dominant_allele_frequency is not 0:
            raise TypeError("dominant_allele_frequency should be float type")
