import threading
from typing import Protocol, Callable, Any
from nicegui import ui


class Signal(Protocol):
	def create_event(self) -> None:
		...

	def set(self, args: dict | None = None) -> None:
		...

	def is_set(self) -> bool:
		...

	def wait(self) -> dict:
		...


class SignalHandler(Protocol):
	def register_signal(self, signal_name: str) -> None:
		...

	def request_signal(self, signal_name: str) -> Signal:
		...

	def send_signal(self, signal_name: str, args: dict | None = None) -> None:
		...

	def clear_signals(self) -> None:
		...


class Controller(Protocol):
	signal_handler: SignalHandler

	def open_page(self, route: str, source: str | None = None,
				  data: dict[str, Any] | None = None) -> None:
		...

	def on_page_open(self, data: dict[str, Any] | None) -> None:
		...

	def start_job(self, source: str, job_fn: Callable[[dict[str, Any]], None], job_data: dict[str, Any]):
		...

	def terminate_job(self) -> None:
		...

	def is_job_running(self) -> bool:
		...

	@staticmethod
	def detect_blastn_installation() -> str | None:
		...

	@staticmethod
	def detect_blastdbcheck_installation() -> str | None:
		...

	@staticmethod
	def check_blast_db(blastdbcheck_binary: str, db_path: str, db_name: str) -> bool:
		...


PageBuilderFn = Callable[[Controller, str, dict[str, Any]], None]
JobFn = Callable[[dict[str, Any], threading.Event, SignalHandler], None]


class PageRouter(Protocol):
	def register_pages(self, pages: dict[str, PageBuilderFn]) -> None:
		...

	def open(self, path: str, source: str | None = None,
			 data: dict[str, Any] | None = None) -> None:
		...

	def frame(self) -> ui.element:
		...

	def build_single_page_frame(self) -> None:
		...
