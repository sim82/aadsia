use aadsia::kdtree::{construct_balanced, KdNode};
use draw::{render, Canvas, Color, Drawing, LineBuilder, Style, SvgRenderer};

fn main() {
    let mut points = [
        [1.0, -1.0],
        [-1.0, 6.0],
        [-0.5, 0.0],
        [2.0, 5.0],
        [0.0, 5.0],
        [2.5, 3.0],
        [-1.0, 1.0],
        [-1.5, -2.0],
        // [1.5, 1.0],
        // [0.0, 0.5],
        // [-2.0, -1.0],
        // [2.0, -3.0],
    ];

    let tree = construct_balanced(&mut points, 0).unwrap();

    let tree = KdNode2::new(
        [0.0, 5.0],
        Some(KdNode2::new(
            [-1.0, 6.0],
            Some(KdNode2::new(
                [-1.0, 1.0],
                Some(KdNode2::new([-1.5, -2.0], None, None, 3)),
                Some(KdNode2::new([-0.5, 0.0], None, None, 3)),
                2,
            )),
            None,
            1,
        )),
        Some(KdNode2::new(
            [1.0, -1.0],
            Some(KdNode2::new(
                [2.0, -5.0],
                None,
                Some(KdNode2::new([2.5, -3.0], None, None, 3)),
                2,
            )),
            None,
            1,
        )),
        0,
    );

    draw_tree(&tree);
}

type KdNode2 = KdNode<2>;

const SCALE: f32 = 50.0;
const XOFFS: f32 = 10.0;
const YOFFS: f32 = 10.0;

fn draw_tree(root: &KdNode2) {
    let mut canvas = Canvas::new(1000, 1000);
    draw_node(root, &mut canvas, -10.0, 10.0, -10.0, 10.0);

    render::save(&canvas, "plot.svg", SvgRenderer::new()).unwrap();
}

fn to_screen_x(x: f32) -> f32 {
    (x + XOFFS) * SCALE
}

fn to_screen_y(y: f32) -> f32 {
    1000.0 - (y + YOFFS) * SCALE
}

fn draw_node(node: &KdNode2, canvas: &mut Canvas, xmin: f32, xmax: f32, ymin: f32, ymax: f32) {
    let x = to_screen_x(node.point[0]);
    let y = to_screen_y(node.point[1]);
    let point = Drawing::new()
        .with_shape(draw::Shape::Circle { radius: 5 })
        .with_xy(x, y)
        .with_style(Style::stroked(1, Color::black()));

    let line_shape = if node.level % 2 == 0 {
        LineBuilder::new(x, to_screen_y(ymin))
            .line_to(x, to_screen_y(ymax))
            .build()
    } else {
        LineBuilder::new(to_screen_x(xmin), y)
            .line_to(to_screen_x(xmax), y)
            .build()
    };

    let line = Drawing::new()
        .with_shape(line_shape)
        .with_style(Style::stroked(1, Color::black()));
    canvas.display_list.add(line);

    canvas.display_list.add(point);

    if let Some(left) = node.left.as_ref() {
        if node.level % 2 == 0 {
            // x line
            draw_node(left, canvas, xmin, node.point[0], ymin, ymax);
        } else {
            // y line
            draw_node(left, canvas, xmin, xmax, ymin, node.point[1]);
        }
    }
    if let Some(right) = node.right.as_ref() {
        if node.level % 2 == 0 {
            // x line
            draw_node(right, canvas, node.point[0], xmax, ymin, ymax);
        } else {
            draw_node(right, canvas, xmin, xmax, node.point[1], ymax);
        }
    }
}
