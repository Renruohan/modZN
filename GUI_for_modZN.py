"""
modZN launcher: lets the user choose between circular and linear assembly and
opens the corresponding graphical interface.

    python GUI_for_modZN.py

On a server without a display use the command-line interfaces directly:

    python modZN_linear.py --help
    python modZN_circular.py --help
"""

import sys

try:
    import FreeSimpleGUI as sg
except ImportError:
    sys.exit("FreeSimpleGUI is not installed (pip install FreeSimpleGUI). "
             "Use the command-line interfaces instead: python modZN_linear.py --help")

import modZN_circular
import modZN_linear


def main() -> None:
    sg.theme("LightBlue")
    sg.theme_button_color(("blue", "#6D9F85"))
    font = ("Helvetica", 14)

    banner = []
    try:
        from Title_png import img as title_image
        banner = [[sg.Image(data=title_image)]]
    except ImportError:
        pass

    layout = banner + [
        [sg.Text("Please choose the way of assembly:", font=font)],
        [sg.Radio("Circular Assembly", "ASSEMBLY", key="-CIRCULAR-", font=font),
         sg.Radio("Linear Assembly", "ASSEMBLY", key="-LINEAR-", font=font)],
        [sg.Button("Primer Design", font=font), sg.Button("Exit", font=font)],
    ]
    window = sg.Window("modZN for primer design", layout, element_justification="center")

    while True:
        event, values = window.read()
        if event in (sg.WIN_CLOSED, "Exit"):
            break
        if event == "Primer Design":
            if values["-CIRCULAR-"]:
                modZN_circular.run_gui()
            elif values["-LINEAR-"]:
                modZN_linear.run_gui()
            else:
                sg.popup_error("Please choose circular or linear assembly.", title="modZN")

    window.close()


if __name__ == "__main__":
    main()
