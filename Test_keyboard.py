# -*- coding: utf-8 -*-
import curses

screen = curses.initscr()
try:
    curses.noecho()
    curses.curs_set(0)
    screen.keypad(1)
    screen.addstr("Press a key")
    event = screen.getch()
finally:
    curses.endwin()

if event == curses.KEY_LEFT or event == curses.KEY_UP:
    print("Left/Up Arrow Key pressed")
elif event == curses.KEY_RIGHT or event == curses.KEY_DOWN:
    print("Right/Down Arrow Key pressed")
else:
    print(event)