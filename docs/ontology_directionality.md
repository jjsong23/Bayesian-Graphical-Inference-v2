# Ontology- and OmniPath-based partial directionality

## Scope

This layer is a conservative first pass for orienting signal propagation. It
does not change the Bayesian probability that an association exists. Instead,
for an undirected edge with probability `p`, it determines whether traversal is
allowed as `A -> B`, `B -> A`, or both. A uniquely disallowed reverse traversal
is stored as zero in the propagation matrix; the allowed traversal retains
exactly `p`.

Two independently auditable sources contribute direction: the ontology-role
rules below and mapped source-target records from the raw mouse OmniPath core
post-translational interaction table. OmniPath direction use has its own GUI
checkbox under the master directionality control and is enabled by default.

The rules are role heuristics, not causal facts entailed by Gene Ontology. In
particular, the kinase/phosphatase-binding rule is an explicit project policy
requested for sensitivity-oriented path inference. It should be revisited when
directed curated interactions, substrate evidence, and biochemical sign are
integrated.

## Active rule catalog

| Source class | Target class | Interpretation |
|---|---|---|
| `ligand` | `receptor` | Ligand acts on receptor. |
| `receptor_regulator` | `receptor` | Receptor regulator acts on receptor. |
| `kinase` | `kinase_phosphatase_binding` | Project policy: kinase-binding protein does not propagate toward the kinase. |
| `phosphatase` | `kinase_phosphatase_binding` | Symmetric project policy for phosphatase-binding proteins. |
| `gtpase_regulator` | `gtpase` | GEF/GAP/regulator acts on GTPase. |
| `cyclase` | `second_messenger` | Cyclase produces cyclic-nucleotide signal. |
| `phospholipase` | `second_messenger` | Phospholipase produces lipid-derived second-messenger signal. |
| `nos` | `second_messenger` | Nitric-oxide synthase produces nitric oxide. |
| `phosphodiesterase` | `second_messenger` | Phosphodiesterase decreases cyclic-nucleotide signal; sign is not yet scored. |
| `second_messenger` | `second_messenger_binding` | Second messenger acts on its binding protein. |

No standalone direction is assigned from `signaling_process`,
`signaling_regulation`, or `adaptor_scaffold`. Those labels establish a role or
proximity, not causal order. Other class combinations are also unresolved
unless they match a listed rule.

## OmniPath direction policy

Every directed OmniPath record whose mouse source and target symbols both map
to the active graph is retained. For each unordered pair:

- if all mapped records point `A -> B`, OmniPath supports only `A -> B`;
- if all mapped records point `B -> A`, OmniPath supports only `B -> A`;
- if records occur in both directions, OmniPath is bidirectional and neither
  traversal is removed.

The OmniPath `consensus_direction` flag, resources, and references are retained
in the audit. The consensus flag is not required for orientation: agreement of
all mapped source-target records for that pair is the operational criterion.
Stimulation and inhibition annotations remain separate from direction and do
not change the edge probability or path score.

## Multi-role and conflict policy

Every class label on both endpoints is considered. A unique ontology
orientation is applied first because these are the established project
constraints. OmniPath then uniquely orients only pairs that ontology left
unresolved. Enabling OmniPath can therefore add disallowed reverse traversals,
but it cannot reopen or reverse an ontology-disallowed traversal.

If OmniPath uniquely agrees with an ontology orientation, that agreement is
recorded. If it uniquely disagrees, the ontology orientation remains in force
and the disagreement is explicitly audited. Opposing ontology rules remain
unresolved unless a unique OmniPath direction resolves them. Bidirectional
OmniPath records do not orient an otherwise unresolved pair. If neither source
supplies direction, the result is `unresolved_no_direction_evidence`; if
direction evidence remains bidirectional, it is
`unresolved_conflicting_directions`. Both traversals remain available in those
unresolved cases.

Every run writes the exact JSON rules, a complete table of all unordered class
pairs, the mapped OmniPath direction table, a pair-level combined
directionality audit, and the partially directed propagation matrix. This makes
the policy independently reviewable and allows rules to be added without
changing Bayesian edge probabilities.
