from typing import Optional
import uuid
import utils


class Controller:

	def __init__(self, model, view) -> None:
		self.model = model
		self.view = view
		self.current_job_id: Optional[str] = None

	def init_new_job(self) -> str:
		job_id = str(uuid.uuid4())
		self.current_job_id = job_id
		return job_id

	def clear_job(self):
		self.current_job_id = None

	def page_allowed(self, page_kwargs: dict) -> bool:
		if not page_kwargs.get('job_id') or not self.current_job_id or \
				page_kwargs.get('job_id') != self.current_job_id:
			return False
		return True

	def handle_fasta_seq_input(self, user_input: str) -> None:
		seq_record = utils.read_fasta_string(user_input)
		self.model.load_sequence_record(seq_record)
		job_id = self.init_new_job()
		print(f'new job: {job_id}')
		self.view.open_page('metainf', job_id=job_id)

	def start_sequence_analysis(self, triplet_settings: dict, regions_settings: dict) -> None:
		self.view.open_page('results', job_id=self.current_job_id)

