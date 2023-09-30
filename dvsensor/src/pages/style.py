from nicegui import ui


class Colors:
	BACKGROUND: str = '#1c1c1c'
	FOREGROUND: str = '#fefefe'
	ACCENT: str = '#d8c5f7'
	RED: str = '#fd6f51'
	GREEN: str = '#019c56'


def set_colors() -> None:
	ui.colors(primary=Colors.BACKGROUND,
			  secondary=Colors.FOREGROUND,
			  accent=Colors.ACCENT)
	ui.query('body').style(f'background-color: {Colors.BACKGROUND}')
