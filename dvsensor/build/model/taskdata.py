from dataclasses import dataclass, field
from typing import Optional
from abc import ABC, abstractmethod
from .rnadata import RNAData


@dataclass
class TaskData(ABC):

	@abstractmethod
	def reload(self) -> None: ...
	
	# @abstractmethod
	# def update(self, key, value) -> None: ...

	# @abstractmethod
	# def query(self, key) -> None: ...


@dataclass
class SingleSeqAnalysisData(TaskData):
	rna_data: RNAData
	len_5UTR: Optional[int] = field(init=False, default=None)
	len_CDS: Optional[int] = field(init=False, default=None)
	len_3UTR: Optional[int] = field(init=False, default=None)
	sensors: list = field(init=False, default_factory=list)

	def __post_init__(self) -> None:
		... # super().__init__()

	def reload(self) -> None:
		original_record = self.rna_data.original_sequence_record
		self.rna_data = RNAData(original_record)
