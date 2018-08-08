import argparse
import csv

def isEmpty(val):
    return val == None or val == "--" or val == "-" or val == "N/A" or val == ""

def checkIfBothFieldsEmpty(field_one, field_two):
    return isEmpty(field_one) and isEmpty(field_two)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sean", help="Sean's priors")
    parser.add_argument("--bx", help="BRCA Exchange priors")
    parser.add_argument("--output", default="priors_concordance.tsv", help="Output file")

    args = parser.parse_args()

    sean_priors = csv.DictReader(open(args.sean, "r"), delimiter="\t")
    bx_priors = csv.DictReader(open(args.bx, "r"), delimiter="\t")

    relevant_fields_from_sean = ["applicable_prior", "protein_prior", "de_novo_prior", "refsplice_prior", "nthgvs", "expert_feature"]
    relevant_fields_from_bx = ["applicablePrior", "proteinPrior", "deNovoDonorPrior", "refDonorPrior", "refAccPrior",
                               "pyhgvs_cDNA", "HGVS_cDNA"]
    ordered_concordance_fieldnames = ["pyhgvs_cDNA", "HGVS_cDNA", "nthgvs", "applicable_prior", "applicablePrior",
                                      "protein_prior", "proteinPrior", "de_novo_prior", "deNovoDonorPrior",
                                      "refsplice_prior", "refDonorPrior", "refAccPrior", "expert_feature", "concordant"]

    concordance = csv.DictWriter(open(args.output, "w"), delimiter="\t",
                             fieldnames=ordered_concordance_fieldnames)

    concordance.writeheader()

    concordance_variants = {}

    # gather all relevant data from bx_priors and store in a dict
    for bx_variant in bx_priors:
        concordance_variant = {}
        for field in relevant_fields_from_bx:
            concordance_variant[field] = bx_variant[field]
        for field in relevant_fields_from_sean:
            concordance_variant[field] = ''

        concordance_variants[bx_variant["pyhgvs_cDNA"].lower().split(":")[1]] = concordance_variant

    # add all relevant data from sean's priors to the dict if variant was found in bx data
    # (combine sean data with bx priors data)
    for sean_variant in sean_priors:
        key = sean_variant["nthgvs"].lower()
        if key in concordance_variants:
            for field in relevant_fields_from_sean:
                concordance_variants[key][field] = sean_variant[field]

    # determine if values are concordant and write to output
    for variant in concordance_variants.values():
        concordant = "True"

        # refsplice_prior corresponds to the union of refDonorPrior and refAccPrior,
        # which in practice are never both set at the same time
        bx_refsplice_prior_equivalent = variant["refAccPrior"] if isEmpty(variant["refDonorPrior"]) else variant["refDonorPrior"]

        # if equivalent fields are equal in value (or both considered empty), priors are concordant
        if ((variant["applicable_prior"] != variant["applicablePrior"] and not checkIfBothFieldsEmpty(variant["applicable_prior"], variant["applicablePrior"]) or
            variant["protein_prior"] != variant["proteinPrior"] and not checkIfBothFieldsEmpty(variant["protein_prior"], variant["proteinPrior"])  or
            variant["de_novo_prior"] != variant["deNovoDonorPrior"] and not checkIfBothFieldsEmpty(variant["de_novo_prior"], variant["deNovoDonorPrior"]) or
            variant["refsplice_prior"] != bx_refsplice_prior_equivalent and not checkIfBothFieldsEmpty(variant["refsplice_prior"], bx_refsplice_prior_equivalent))):
            concordant = "False"

        variant["concordant"] = concordant

        # only output variants that have data from sean
        if not isEmpty(variant['nthgvs']):
            concordance.writerow(variant)


if __name__ == "__main__":
    main()
