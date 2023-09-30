import logging
import threading
from dataclasses import dataclass, field


@dataclass
class Signal:
	name: str
	requested: bool = field(default=False)
	event: threading.Event | None = field(default=None, init=False)
	args: dict = field(default_factory=dict, init=False)
	_event_created: threading.Event = field(default_factory=threading.Event, init=False)

	def __post_init__(self) -> None:
		if not self.requested:
			self.create_event()

	def create_event(self) -> None:
		if not self.event:
			self.event = threading.Event()
			self._event_created.set()

	def set(self, args: dict | None = None) -> None:
		if self.event:
			if args:
				self.args = args
			self.event.set()

	def is_set(self) -> bool:
		if self.event:
			return self.event.is_set()
		return False

	def wait(self) -> dict:
		if not self.event:
			self._event_created.wait()
		self.event.wait()
		return self.args


class SignalHandler:

	def __init__(self) -> None:
		self._signals: dict[str, Signal] = {}

	def register_signal(self, signal_name: str) -> None:
		if signal_name in self._signals:
			if self._signals[signal_name].requested:
				self._signals[signal_name].create_event()
				logging.debug(f"requested signal '{signal_name}' registered")
			else:
				logging.warning(f"could not register signal {signal_name}: already registered")
		else:
			self._signals[signal_name] = Signal(signal_name)
			logging.debug(f"signal '{signal_name}' registered")

	def request_signal(self, signal_name: str) -> Signal:
		if signal_name in self._signals:
			return self._signals[signal_name]
		signal = Signal(signal_name, requested=True)
		self._signals[signal_name] = signal
		return signal

	def send_signal(self, signal_name: str, args: dict | None = None) -> None:
		if signal_name in self._signals:
			self._signals[signal_name].set(args)
			logging.debug(f"signal '{signal_name}' sent")
		else:
			logging.warning(f"could not send signal '{signal_name}': not registered")

	def clear_signals(self) -> None:
		self._signals = {}
