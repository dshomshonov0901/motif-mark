import cairo
import os

WIDTH, HEIGHT = 500, 400

surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
ctx = cairo.Context(surface)

# white background
ctx.set_source_rgb(1, 1, 1)
ctx.paint()

# black line
ctx.set_source_rgb(0, 0, 0)
ctx.set_line_width(10)
ctx.move_to(450, 50)
ctx.line_to(350, 350)
ctx.stroke()

# rectangle outline
ctx.set_line_width(6)
ctx.rectangle(100, 100, 200, 120)
ctx.stroke()

filename = "pycairo_basics.png"
surface.write_to_png(filename)

