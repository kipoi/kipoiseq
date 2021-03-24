import abc
from typing import List, Union

import pandas as pd
import numpy as np

from kipoiseq import Interval

import logging

from kipoiseq.extractors.multi_interval import BaseMultiIntervalFetcher

log = logging.getLogger(__name__)

__all__ = [
    "GTFMultiIntervalFetcher",
    "CDSFetcher",
    "UTRFetcher",
    "gtf_row2interval"
]


def _filter_valid_transcripts(gtf_df):
    if 'transcript_support_level' in gtf_df:
        gtf_df = gtf_df[~gtf_df.transcript_support_level.isnull()]
        gtf_df = gtf_df[gtf_df.transcript_support_level != 'NA']
        gtf_df.transcript_support_level = gtf_df.transcript_support_level.astype(
            int)
        gtf_df = gtf_df[~gtf_df.transcript_support_level.isna()]
    else:
        raise ValueError('Transcript support level not in gtf. '
                         'Skipping the associated filters.')
    return gtf_df


def _get_biotype_str(df):
    if 'transcript_biotype' in df:
        return 'transcript_biotype'
    elif 'gene_biotype' in df:
        return 'gene_biotype'
    elif 'gene_type' in df:
        return 'gene_type'
    else:
        raise ValueError('Cannot obtain `biotype_str` from gtf file')


def _filter_biotype(gtf_df, value):
    """
    Gets the biotype column and checks whether it equals `value`.
    :param gtf_df: genome annotation GTF dataframe
    :param value: The value to check for
    :return: filtered dataframe
    """
    biotype_str = _get_biotype_str(gtf_df)
    return gtf_df.query(
        "{key} == '{value}'".format(key=biotype_str, value=value)
    )


def _filter_biotype_proteincoding(gtf_df):
    return _filter_biotype(gtf_df, value="protein_coding")


def _filter_tag(gtf_df, regex_contains):
    return gtf_df.query(
        "tag.str.contains('{}').fillna(False).values".format(regex_contains)
    )


def gtf_row2interval(row, interval_attrs: Union[str, List[str]] = None):
    """
    Convert gtf row object into interval class.
    """
    if interval_attrs:
        # check type
        if not isinstance(interval_attrs, list):
            interval_attrs = [interval_attrs]

        interval_attrs = {i: row[i] for i in interval_attrs if i in row}

    return Interval(str(row.Chromosome),
                    int(row.Start),
                    int(row.End),
                    strand=str(row.Strand),
                    attrs=interval_attrs)


class GTFMultiIntervalFetcher(BaseMultiIntervalFetcher):
    df: pd.DataFrame

    def __init__(self, region_df: pd.DataFrame, keep_attrs=None, on_error_warn=True, **kwargs):
        """
        Uses a a pyranges/GTF-like dataframe as source for intervals.

        :param region_df:
            GTF dataframe as given by PyRanges containing at least the following columns:

                - Chromosome
                - Start
                - End
                - Strand

            Will use `region_df.index` as source for all keys.
            By assigning multiple rows with the same index, these rows will be converted into a list of intervals when
            requesting a key.
        :param keep_attrs: (list of) attributes to keep during conversion of dataframe columns to intervals
        :param on_error_warn: should errors raise an exception or just log a warning?
        """
        self.df = region_df
        self.keep_attrs = keep_attrs
        self.on_error_warn = on_error_warn

    def keys(self) -> List[object]:
        return self.df.index.unique()

    def _gtf_df2interval(self, intervals_df: pd.DataFrame) -> List[Interval]:
        """
        Convert a pyranges/GTF-like dataframe into a list of intervals

        :param intervals_df: GTF dataframe as given by PyRanges containing at least the following columns:
            - Chromosome
            - Start
            - End
            - Strand
        :return: list of intervals
        """
        return [
            gtf_row2interval(row, self.keep_attrs) for i, row in intervals_df.sort_values("Start").iterrows()
        ]

    def get_intervals(self, key) -> List[Interval]:
        """
        Returns a list of Intervals on request of some key
        :param key: The identifier for the requested intervals
        :return: list of intervals
        """
        # make sure that key is a list;
        # otherwise, the resulting subset could be a series instead of a dataframe
        if not isinstance(key, list):
            key = [key]

        intervals_df = self.df.loc[key]

        check_strand = np.size(np.unique(intervals_df.Strand)) == 1
        if not check_strand:
            if self.on_error_warn:
                log.warning("Error while processing item '%s': strands differ", key)
            else:
                raise ValueError("Error while processing item '%s': strands differ" % key)

        return self._gtf_df2interval(intervals_df)


class CDSFetcher(GTFMultiIntervalFetcher):

    def __init__(
            self,
            gtf_file,
            filter_valid_transcripts=True,
            filter_biotype=True,
            filter_tag=True,
            on_error_warn=True,
    ):
        """
        Interval fetcher for coding sequences
        :param gtf_file: input GTF file
        """
        # TODO: add docs for parameters
        self.gtf_file = str(gtf_file)

        df = self._read_cds(
            self.gtf_file,
            filter_valid_transcripts=filter_valid_transcripts,
            filter_biotype=filter_biotype,
            filter_tag=filter_tag,
            on_error_warn=on_error_warn
        )

        super().__init__(
            region_df=df,
            keep_attrs="tag",
            on_error_warn=on_error_warn
        )

    @property
    def transcripts(self):
        return self.keys()

    @staticmethod
    def _read_cds(
            gtf_file,
            filter_valid_transcripts=False,
            filter_biotype=False,
            filter_tag=False,
            duplicate_attr=None,
            on_error_warn=True,
    ):
        """
        Read, extract and filter valid cds from the given gtf_file
        :param gtf_file: path to the GTF file
        """
        import pyranges

        if duplicate_attr == None:
            # One row in the GTF file can have multiple tags;
            # therefore, to filter them we have to allow duplicate attrs.
            duplicate_attr = filter_tag

        df = pyranges.read_gtf(gtf_file, as_df=True, duplicate_attr=duplicate_attr)

        cds = CDSFetcher.get_cds_from_gtf(
            df,
            filter_valid_transcripts=filter_valid_transcripts,
            filter_biotype=filter_biotype,
            filter_tag=filter_tag,
            on_error_warn=on_error_warn
        )

        cds = cds.set_index("transcript_id")
        return cds

    @staticmethod
    def get_cds_from_gtf(
            df,
            filter_valid_transcripts=False,
            filter_biotype=False,
            filter_tag=False,
            on_error_warn=True,
    ):
        """
        Create DataFrame with valid cds

        :param df: the GTF dataframe
        :param filter_valid_transcripts: Filter for "transcript_support_level" column in GTF
        :param filter_biotype: Filter for "[gene|transcript]_biotype == 'protein_coding'" in GTF
        :param filter_tag: Filter for 'basic' or 'CCDS' tag
        """
        cds = df.query("(Feature == 'CDS') | (Feature == 'CCDS')")

        if filter_biotype:
            try:
                cds = _filter_biotype_proteincoding(cds)
            except ValueError as e:
                if on_error_warn:
                    log.warning("Error during filtering biotype: %s", e)
                else:
                    raise e
        if filter_valid_transcripts:
            try:
                cds = _filter_valid_transcripts(cds)
            except ValueError as e:
                if on_error_warn:
                    log.warning("Error during filtering for valid transcripts: %s", e)
                else:
                    raise e
        if filter_tag:
            try:
                cds = _filter_tag(cds, regex_contains='basic|CCDS')
            except ValueError as e:
                if on_error_warn:
                    log.warning("Error during filtering for tag: %s", e)
                else:
                    raise e

        return cds


class UTRFetcher(GTFMultiIntervalFetcher):

    def __init__(
            self,
            gtf_file,
            feature_type="5UTR",
            infer_from_cds=False,
            on_error_warn=True,
    ):
        """
        :param gtf_file: path to the GTF file
        :param feature_type: type of the feature that will be filtered for. In general '5UTR' or '3UTR'.
        :param infer_from_cds: Substract the CDS from the exon regions to infer the UTR regions.
            Will use 'feature_type' to decide whether '5UTR' or '3UTR' should be returned.
        :param on_error_warn: Do not break on error; instead throw warning.
        """
        self.gtf_file = str(gtf_file)

        df = self._read_utr(
            self.gtf_file,
            feature_type=feature_type,
            infer_from_cds=infer_from_cds,
            on_error_warn=on_error_warn
        )

        super().__init__(
            region_df=df,
            keep_attrs=None,
            on_error_warn=on_error_warn
        )

    @property
    def transcripts(self):
        return self.keys()

    @staticmethod
    def _read_utr(
            gtf_file,
            feature_type="5UTR",
            infer_from_cds=False,
            on_error_warn=True,
    ) -> pd.DataFrame:
        """
        Read, extract and filter valid UTRs from the given gtf_file
        :param gtf_file: path to the GTF file
        :param feature_type: type of the feature that will be filtered for. In general '5UTR' or '3UTR'.
        :param infer_from_cds: Substract the CDS from the exon regions to infer the UTR regions.
            Will use 'feature_type' to decide whether '5UTR' or '3UTR' should be returned.
        :param on_error_warn: Do not break on error; instead throw warning.
        """
        import pyranges

        df = pyranges.read_gtf(gtf_file, as_df=True)

        utr_df = UTRFetcher.get_utr_from_gtf(
            df,
            feature_type=feature_type,
            infer_from_cds=infer_from_cds,
            on_error_warn=on_error_warn
        )

        utr_df = utr_df.set_index("transcript_id")
        return utr_df

    @staticmethod
    def get_utr_from_gtf(
            df,
            feature_type="5UTR",
            infer_from_cds=False,
            on_error_warn=True,
    ) -> pd.DataFrame:
        """
        Create DataFrame with valid UTRs
        :param df: the GTF dataframe
        :param feature_type: type of the feature that will be filtered for. In general '5UTR' or '3UTR'.
        :param infer_from_cds: Substract the CDS from the exon regions to infer the UTR regions.
            Will use 'feature_type' to decide whether '5UTR' or '3UTR' should be returned.
        :param on_error_warn: Do not break on error; instead throw warning.
        """
        if infer_from_cds:
            # get start and end of cds for each transcript
            cds = (
                CDSFetcher.get_cds_from_gtf(df=df, on_error_warn=on_error_warn)
                .groupby('transcript_id')
                .agg({'Start': min, 'End': max})
            )

            # join cds start and end to utr df
            utr_df = (
                df.query("Feature == 'transcript'")
                .set_index('transcript_id')
                .join(cds, rsuffix="_cds")
                .dropna(subset=['Start_cds', 'End_cds'], axis=0)
            )

            if feature_type.upper() == "5UTR":
                utr_df['Start'] = np.where(utr_df['Strand'] == '+', utr_df['Start'].astype("int"), utr_df['End_cds'].astype("int"))
                utr_df['End'] = np.where(utr_df['Strand'] == '+', utr_df['Start_cds'].astype("int"), utr_df['End'].astype("int"))
                utr_df['Feature'] = pd.Categorical("5UTR", categories = utr_df['Feature'])
            if feature_type.upper() == "3UTR":
                utr_df['Start'] = np.where(utr_df['Strand'] == '+', utr_df['End_cds'].astype("int"), utr_df['Start'].astype("int"))
                utr_df['End'] = np.where(utr_df['Strand'] == '+', utr_df['End'].astype("int"), utr_df['Start_cds'].astype("int"))
                utr_df['Feature'] = pd.Categorical("3UTR", categories = utr_df['Feature'])

            utr_df.drop(['Start_cds', 'End_cds'], axis=1, inplace=True)
            return utr_df.reset_index()
        else:
            utr_df = df.query("Feature == '{feature_type}'".format(feature_type=feature_type))

        return utr_df
