from dataclasses import dataclass, field
from typing import Optional
from abc import ABC, abstractmethod
from .rnadata import RNAData
import asyncio


@dataclass
class TaskData(ABC):

	@abstractmethod
	def reload(self) -> None: ...
	
	@abstractmethod
	def update(self, key, value) -> None: ...

	@abstractmethod
	async def query(self, key): ...


@dataclass
class SingleSeqAnalysisData(TaskData):
	rna_data: RNAData
	len_5UTR: Optional[int] = field(init=False)
	len_CDS: Optional[int] = field(init=False)
	len_3UTR: Optional[int] = field(init=False)
	sensors: list = field(init=False, default_factory=list)
	initialized: dict = field(init=False)

	def __post_init__(self) -> None:
		self.initialized = {}

	def reload(self) -> None:
		original_record = self.rna_data.original_sequence_record
		self.rna_data = RNAData(original_record)

	def update(self, key, value) -> None:
		if not hasattr(self, key) and key in self.initialized.keys():
			self.initialized[key].set()
			del self.initialized[key]
		setattr(self, key, value)

	async def query(self, key):

		if hasattr(self, key):
			return getattr(self, key)
		initialized = asyncio.Event()
		self.initialized[key] = initialized
		await initialized.wait()

		return getattr(self, key)



