from dataclasses import dataclass, field

@dataclass
class ProteinIdentity:
    # reprisents each protien
    accession: str
    name: str | None = None

@dataclass
class StructureSummary:
    # reprisents what is physically present in the PDB file fetched by ppinsight 
    pdb_path: str | None = None
    chains: list[str] = field(default_factory=list)
    residue_numbers: dict[str, list[int]] = field(default_factory=dict)
    hetero_residues: dict[str, list[int]] = field(default_factory=dict)

@dataclass
class ComplexEvidence:
    pdb_id: str
    source: str
    chain_a: str | None = None
    chain_b: str | None = None
    pair_relevant: bool = False

@dataclass
class EvidenceRecord:
    source: str
    evidence_type: str
    description: str
    protein: str | None = None
    residues: list[int] = field(default_factory=list)
    region_start: int | None = None
    region_end: int | None = None

    pair_relevant: bool = False # false default, true if the evidence is relevant to the pair of proteins in question
    source_id: str | None = None

@dataclass
class PairEvidence:
    protein_a: ProteinIdentity
    protein_b: ProteinIdentity
    # known region: state 1
    interface_regions_a: list[EvidenceRecord] = field(default_factory=list)
    interface_regions_b: list[EvidenceRecord] = field(default_factory=list)
    # known residues: state 2
    interface_residues_a: list[EvidenceRecord] = field(default_factory=list)
    interface_residues_b: list[EvidenceRecord] = field(default_factory=list)
    # known complex: state 3
    known_complex: ComplexEvidence | None = None
    # starting pose: state 3
    starting_pose: str | None = None
    evidence_records: list[EvidenceRecord] = field(default_factory=list)

@dataclass
class KnowledgeState:
    level: int
    name: str
    recommended_future_strategy: str
    reason: str

@dataclass
class EvidenceReport:
    pair_evidence: PairEvidence
    knowledge_state: KnowledgeState