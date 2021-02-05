# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import collections

from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import QIIME1DemuxDirFmt
import biom
import skbio.io

from ._cluster_features import run_command


def _parse_uc(fh):
    """ This function is copied from the biom-format project temporarily to
        fix a bug in handling of sample ids with underscores in them
        (https://github.com/biocore/biom-format/issues/758). This code
        will be contributed back upstream to the biom-format project, and
        will be removed from this plugin when a biom-format release is
        available that contains this fix.

       Create a Table object from a uclust/usearch/vsearch uc file.
        Parameters
        ----------
        fh : file handle
            The ``.uc`` file to be parsed.
        Returns
        -------
        biom.Table : The resulting BIOM table.
        Raises
        ------
        ValueError
            If a sequence identifier is encountered that doesn't have at least
            one underscore in it (see Notes).
        Notes
        -----
        This function assumes sequence identifiers in this file are formated as
        ``<sample-id>_<sequence-id>``. Everything before the last underscore
        will be used as the sample identifier in the resulting ``Table``.
        The information after the last underscore is not used directly, though
        the full identifiers of seeds will be used as the observation
        identifier in the resulting ``Table``.
    """
    data = collections.defaultdict(int)
    sample_idxs = {}
    sample_ids = []
    observation_idxs = {}
    observation_ids = []
    # The types of hit lines we need here are hit (H), seed (S) and
    # library seed (L). Store these in a set for quick reference.
    line_types = set('HSL')
    for line in fh:
        # determine if the current line is one that we need
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')

        line_type = fields[0]
        if line_type not in line_types:
            continue

        # grab the fields we care about
        observation_id = fields[9].split()[0]
        query_id = fields[8].split()[0]

        if observation_id == '*':
            # S and L lines don't have a separate observation id
            observation_id = query_id

        # get the index of the current observation id, or create it if it's
        # the first time we're seeing this id
        if observation_id in observation_idxs:
            observation_idx = observation_idxs[observation_id]
        else:
            observation_idx = len(observation_ids)
            observation_ids.append(observation_id)
            observation_idxs[observation_id] = observation_idx

        if line_type == 'H' or line_type == 'S':
            # get the sample id
            try:
                # the following line is modified from biom-format 2.1.6 to
                # find the last underscore rather than the first
                underscore_index = query_id.rindex('_')
            except ValueError:
                raise ValueError(
                 "A query sequence was encountered that does not have an "
                 "underscore. An underscore is required in all query "
                 "sequence identifiers to indicate the sample identifier.")

            # get the sample id and its index, creating the index if it is the
            # first time we're seeing this id
            sample_id = query_id[:underscore_index]
            if sample_id in sample_idxs:
                sample_idx = sample_idxs[sample_id]
            else:
                sample_idx = len(sample_ids)
                sample_idxs[sample_id] = sample_idx
                sample_ids.append(sample_id)
            # increment the count of the current observation in the current
            # sample by one.
            data[(observation_idx, sample_idx)] += 1
        else:
            # nothing else needs to be done for 'L' records
            pass
    return biom.Table(data, observation_ids=observation_ids,
                      sample_ids=sample_ids)


def dereplicate_sequences(sequences: QIIME1DemuxDirFmt,
                          derep_prefix: bool = False
                          ) -> (biom.Table, DNAFASTAFormat):
    dereplicated_sequences = DNAFASTAFormat()
    with tempfile.NamedTemporaryFile(mode='w+') as out_uc:
        seqs_fp = '%s/seqs.fna' % str(sequences)
        cmd = ['vsearch',
               '--derep_fulllength', seqs_fp,
               '--output', str(dereplicated_sequences),
               '--relabel_sha1', '--relabel_keep',
               '--uc', out_uc.name,
               '--qmask', 'none',
               '--xsize',
               '--minseqlength', '1',
               '--fasta_width', '0']
        if derep_prefix:
            cmd[1] = '--derep_prefix'
        run_command(cmd)
        out_uc.seek(0)
        table = _parse_uc(out_uc)
    id_map = {e.metadata['description']: e.metadata['id']
              for e in skbio.io.read(str(dereplicated_sequences),
                                     constructor=skbio.DNA,
                                     format='fasta')}
    table.update_ids(id_map=id_map, axis='observation')
    return table, dereplicated_sequences
