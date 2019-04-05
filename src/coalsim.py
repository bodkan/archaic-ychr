#!/usr/bin/env python3

import argparse

import pandas as pd
import msprime as msp


GEN_TIME = 25
Y_Ne = 0.5 * 0.25


def years_to_gen(y):
    return int(y  / GEN_TIME)


def define_samples(pop_params):
    """Generate list of sample definitions for msprime."""
    sample_names = []
    for i, pop in enumerate(pop_params):
        times = [years_to_gen(t) for t in pop_params[pop]["t_sample"]]
        sample_names.extend([msp.Sample(population=i, time=t) for t in times])
    return sample_names


def all_inds(pop_params):
    """Generate list of all simulated sample names."""
    sample_names = []
    for p in pop_params:
        n_pop = len(pop_params[p]["t_sample"])
        sample_names.extend([f"{p}{i}" for i in range(n_pop)])
    return sample_names


parser = argparse.ArgumentParser()

parser.add_argument("--Ne-nonafr", type=int, default=5000)
parser.add_argument("--Ne-afr", type=int, default=10000)
parser.add_argument("--Ne-arch", type=int, default=1000)

parser.add_argument("--seq-len", help="Sequence length", type=int, default=7_000_000)
parser.add_argument("--mut-rate", help="Mutation rate [per bp per year]", type=float, default=7.6e-10)

parser.add_argument("--ui-age", nargs="+", type=int, help="Age of the Ust-Ishim individual [years BP]")
parser.add_argument("--arch-ages", nargs="+", type=int, help="Ages of archaic samples [years BP]", required=True)
parser.add_argument("--neur", type=int, help="Number of European chromosomes to simulate", required=True)
parser.add_argument("--nasn", type=int, default=0, help="Number of Asian chromosomes to simulate")
parser.add_argument("--nafr", type=int, default=1, help="Number of African chromosomes to simulate")

parser.add_argument("--split-chimp", type=int, default=6_000_000)
parser.add_argument("--split-arch", type=int, required=True)
parser.add_argument("--split-afr", type=int, required=True)
parser.add_argument("--split-asn", type=int, default=45_000)

parser.add_argument("--scale", type=float, default=1.0)

parser.add_argument("--output", metavar="FILE", help="Output filename")
parser.add_argument("--format", choices=["snp", "fa"], required=True)

parser.add_argument("--debug", action="store_true", help="Print debugging info", default=False)

args = parser.parse_args()
# args = parser.parse_args("--debug --arch-ages 13000 --split-arch 650_000 --split-afr 250_000 "
# "--ui-age 45000 --neur 5 --nafr 5 --nasn 5 --output out --format snp".split())

pop_params = {
    "chimp": {"id": 0, "Ne": Y_Ne * 10000, "t_sample": 1 * [0], "t_split" : args.split_chimp},
    "arch": {"id": 1, "Ne": Y_Ne * args.Ne_arch,  "t_sample": args.arch_ages, "t_split": args.split_arch},
    "afr": {"id": 2, "Ne": Y_Ne * args.Ne_afr, "t_sample": args.nafr * [0]},
    "eur": {"id": 3, "Ne": Y_Ne * args.Ne_nonafr, "t_sample": args.ui_age + args.neur * [0], "t_split": args.split_afr},
    "asn": {"id": 4, "Ne": Y_Ne * args.Ne_nonafr, "t_sample": args.nasn * [0], "t_split": args.split_asn}
}

id_chimp, id_arch, id_afr, id_eur, id_asn = [pop_params[p]["id"] for p in pop_params]

t_chimp, t_arch, t_eur, t_asn = [
    years_to_gen(pop_params[p]["t_split"])
    for p in ["chimp", "arch", "eur", "asn"]
]

demography = [
    # EUR-ASN split
    msp.MassMigration(t_asn, id_asn, id_eur, 1.0),

    # population size during the bottleneck
    msp.PopulationParametersChange(
        time=years_to_gen(60000),
        initial_size=Y_Ne * 2000,
        population_id=id_eur
    ),

    # population size during the bottleneck
    msp.PopulationParametersChange(
        time=years_to_gen(80000),
        initial_size=Y_Ne * args.Ne_afr,
        population_id=id_eur
    ),

    # out of Africa migration
    msp.MassMigration(t_eur, id_eur, id_afr, 1.0),

    # split of the archaic lineage from modern humans
    msp.MassMigration(t_arch, id_arch, id_afr, 1.0),

    # split of chimps from the human lineage
    msp.MassMigration(t_chimp, id_chimp, id_afr, 1.0),
]

demography.sort(key=lambda event: event.time)

pop_config = [msp.PopulationConfiguration(initial_size=pop_params[p]["Ne"]) for p in pop_params]
samples = define_samples(pop_params)

if args.debug:
    print("Population identifiers:")
    for p, i in pop_params.items():
        print(f"    {p}: {i['id']}")
    msp.DemographyDebugger(
        population_configurations=pop_config,
        demographic_events=demography
    ).print_history()

ts = msp.simulate(
    length=args.seq_len,
    mutation_rate=args.mut_rate * GEN_TIME,
    recombination_rate=0,
    samples=samples,
    population_configurations=pop_config,
    demographic_events=demography
)


snps = pd.DataFrame(
    ts.genotype_matrix(),
    columns=all_inds(pop_params),
    index=(v.position for v in ts.variants())
).rename(columns={"eur0": f"ustishim"})

if args.format == "snp":
  snps.to_csv(args.output + ".snp.gz", sep="\t", index_label="pos")
else:
  with open(args.output + ".fa", "w") as fa:
      for s in snps.columns:
          print(f">{s}", file=fa)
          print("".join("A" if i == 0 else "T" for i in snps[s]), file=fa)
