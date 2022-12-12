import heatrapy as htp

# Create a 2D model
example = htp.SingleObject2D(
    293,
    material="Cu",
    boundaries=(300, 0, 0, 0),
    size=(20, 20)
)

example.change_material(
     material='Al',
     shape='square',
     initial_point=(5, 2),
     length=(4, 4)
)

example.change_material(
        material='water',
        shape='circle',
        initial_point=(10, 10),
        length=4
)

example.compute(100, 10)

example.show_figure("temperatre")
example.show_figure("materials")
example.show_figure("state")
example.show_figure("Q")
example.show_figure("Q0")