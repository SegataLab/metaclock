import re

class ContigsType(object):
    def __init__(self, reference_genome, samples, **kwargs):
        self.input_type = 'contigs'
        self.reference_genome = reference_genome
        self.samples = samples
        self.intermediate = kwargs.get('intermediate', '')

        if 'parameter_set' not in kwargs:
            raise TypeError('parameter_set is not provided')
        parameter_set = kwargs['parameter_set']

        self.param_set = {}
        self.param_set['homolog_length'] = parameter_set.get(
                                            'homolog_length', 500)
        self.param_set['homolog_identity'] = parameter_set.get(
                                            'homolog_identity', 95.0)
        self.param_set['blastn_threads'] = parameter_set.get('blastn_threads', 6)

        self.validate()

    def validate(self):
        if type(self.reference_genome) is not str:
            raise TypeError("reference_genome should be str type")

        # if type(self.samples) is not str:
        #     raise TypeError("samples should be str type")

        if type(self.param_set['homolog_length']) is not int:
            raise TypeError("homolog_length should be an integer")

        if type(self.param_set['homolog_identity']) is not float:
            raise TypeError("homolog_identity should be an float number")

        if type(self.param_set['blastn_threads']) is not int:
            raise TypeError("blastn_threads should be an integer")


class ReadsType(object):
    def __init__(self, reference_genome, samples):
        self.input_type = 'reads'
        self.reference_genome = reference_genome
        self.samples = samples

        # self.validate()

    def validate(self):
        if type(self.reference_genome) is not str:
            raise TypeError("reference_genome should be str type")


class AncientReadsType(ReadsType):
    def __init__(self, **kwargs):
        super().__init__(kwargs['reference_genome'], kwargs['samples'])

        self.age_type = 1
        self.intermediate = kwargs.get('intermediate', 'intermediates')

        if 'parameter_set' not in kwargs:
            raise TypeError('parameter_set is not provided')
        parameter_set = kwargs['parameter_set']

        self.param_set = {}
        self.param_set['search_report_mode'] = parameter_set.get('search_report_mode', '-k,5')
        self.param_set['bowtie2_threads'] = parameter_set.get('bowtie2_threads', 1)
        self.param_set['minimum_mapping_quality'] = parameter_set.get('minimum_mapping_quality', 30)
        self.param_set['minimum_mapping_length'] = parameter_set.get('minimum_mapping_length', 30)
        self.param_set['maximum_snp_edit_distance'] = parameter_set.get('maximum_snp_edit_distance', 0.03)
        self.param_set['nproc'] = parameter_set.get('nproc', 1)
        self.param_set['minimum_coverage'] = parameter_set.get('minimum_coverage', 5)
        self.param_set['trim_distance'] = parameter_set.get('trim_distance', '5:5')
        self.param_set['dominant_allele_frequency'] = parameter_set.get('dominant_allele_frequency', 0.8)
        self.param_set['output_trimmed_reads'] = parameter_set.get('output_trimmed_reads', 0)

        self.validate()

    def validate(self):
        super().validate()

        if type(self.intermediate) is not str:
            raise TypeError("intermediate should be str type")

        if type(self.param_set['search_report_mode']) is not str:
            raise TypeError("search_report_mode should be str type")

        if type(self.param_set['bowtie2_threads']) is not int:
            raise TypeError("bowtie2_threads should be int type")

        if type(self.param_set['minimum_mapping_quality']) is not int:
            raise TypeError("minimum_mapping_quality should be int type")

        if type(self.param_set['minimum_mapping_length']) is not int:
            raise TypeError("minimum_mapping_length should be int type")

        if type(self.param_set['maximum_snp_edit_distance']) is not float:
            raise TypeError("maximum_snp_edit_distance should be float type")

        if type(self.param_set['nproc']) is not int:
            raise TypeError("nproc should be int type")

        if type(self.param_set['minimum_coverage']) is not int:
            raise TypeError("minimum_coverage should be int type")

        if type(self.param_set['trim_distance']) is not str:
            raise TypeError("trim_distance should be str type")

        if not re.match('[0-9]:[0-9]', self.param_set['trim_distance']):
            raise TypeError("trim_distance should ' \
                            'match re pattern `[0-9]:[0-9]`")

        if type(self.param_set['dominant_allele_frequency']) is not float \
            and self.param_set['dominant_allele_frequency'] is not 0:
            raise TypeError("dominant_allele_frequency should be float type")

        if type(self.param_set['output_trimmed_reads']) is not int:
            raise TypeError("output_trimmed_reads should be int type")


class ModernReadsType(ReadsType):
    def __init__(self, **kwargs):
        super().__init__(kwargs['reference_genome'], kwargs['samples'])

        self.age_type = 2
        self.intermediate = kwargs.get('intermediate', 'intermediates')

        if 'parameter_set' not in kwargs:
            raise TypeError('parameter_set is not provided')
        parameter_set = kwargs['parameter_set']

        self.param_set = {}
        self.param_set['search_report_mode'] = parameter_set.get('search_report_mode', '-k,5')
        self.param_set['bowtie2_threads'] = parameter_set.get('bowtie2_threads', 1)
        self.param_set['minimum_mapping_quality'] = parameter_set.get('minimum_mapping_quality', 30)
        self.param_set['minimum_mapping_length'] = parameter_set.get('minimum_mapping_length', 30)
        self.param_set['maximum_snp_edit_distance'] = parameter_set.get('maximum_snp_edit_distance', 0.03)
        self.param_set['nproc'] = parameter_set.get('nproc', 1)
        self.param_set['minimum_coverage'] = parameter_set.get('minimum_coverage', 5)
        self.param_set['dominant_allele_frequency'] = parameter_set.get('dominant_allele_frequency', 0.8)

        self.validate()

    def validate(self):
        super().validate()

        if type(self.intermediate) is not str:
            raise TypeError("intermediate should be str type")

        if type(self.param_set['search_report_mode']) is not str:
            raise TypeError("search_report_mode should be str type")

        if type(self.param_set['bowtie2_threads']) is not int:
            raise TypeError("bowtie2_threads should be int type")

        if type(self.param_set['minimum_mapping_quality']) is not int:
            raise TypeError("minimum_mapping_quality should be int type")

        if type(self.param_set['minimum_mapping_length']) is not int:
            raise TypeError("minimum_mapping_length should be int type")

        if type(self.param_set['maximum_snp_edit_distance']) is not float:
            raise TypeError("maximum_snp_edit_distance should be float type")

        if type(self.param_set['nproc']) is not int:
            raise TypeError("nproc should be int type")

        if type(self.param_set['minimum_coverage']) is not int:
            raise TypeError("minimum_coverage should be int type")

        if type(self.param_set['dominant_allele_frequency']) is not float \
            and self.param_set['dominant_allele_frequency'] is not 0:
            raise TypeError("dominant_allele_frequency should be float type")

