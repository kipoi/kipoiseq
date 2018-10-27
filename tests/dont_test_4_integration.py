import kipoi
import kipoi_veff
import kipoi_veff.snv_predict as sp
import pytest
from kipoi_veff.utils.generic import ModelInfoExtractor
from kipoiseq.dataloaders.sequence import IntervalSeqDl
import os
from kipoi.pipeline import install_model_requirements

INSTALL_REQ = False


def test_deepsea():
    model = kipoi.get_model("DeepSEA/variantEffects")
    mie = ModelInfoExtractor(model, IntervalSeqDl)


def test_var_eff_pred_varseq(tmpdir):
    model_name = "DeepSEA/variantEffects"
    if INSTALL_REQ:
        install_model_requirements(model_name, "kipoi", and_dataloaders=True)
    #
    model = kipoi.get_model(model_name, source="kipoi")
    # The preprocessor
    Dataloader = IntervalSeqDl
    #
    dataloader_arguments = {"intervals_file": "example_files/intervals.bed",
                            "fasta_file": "example_files/hg38_chr22.fa",
                            "required_seq_len": 1000, "alphabet_axis": 1, "dummy_axis": 2, "label_dtype": str}
    dataloader_arguments = {k: model.source_dir + "/" + v if isinstance(v, str) else v for k, v in
                            dataloader_arguments.items()}

    vcf_path = "tests/data/variants.vcf"
    out_vcf_fpath = str(tmpdir.mkdir("variants_generated", ).join("out.vcf"))
    #
    vcf_path = kipoi_veff.ensure_tabixed_vcf(vcf_path)
    model_info = kipoi_veff.ModelInfoExtractor(model, Dataloader)
    writer = kipoi_veff.VcfWriter(model, vcf_path, out_vcf_fpath, standardise_var_id=True)
    vcf_to_region = kipoi_veff.SnvCenteredRg(model_info)
    res = sp.predict_snvs(model, Dataloader, vcf_path, dataloader_args=dataloader_arguments,
                          batch_size=32,
                          vcf_to_region=vcf_to_region,
                          sync_pred_writer=writer)
    writer.close()
    assert os.path.exists(out_vcf_fpath)
