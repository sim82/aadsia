use aadsia::sstree::{SsNode, SsTree};
use draw::{render, Canvas, Color, Drawing, LineBuilder, Shape, Style, SvgRenderer, RGB};
use rand::{random, Rng};

fn main() {
    let mut tree = SsTree::<2, 8>::new(4);
    let mut rng = rand::thread_rng();

    let mut points = Vec::new();
    for _ in 0..12 {
        let point = [rng.gen_range(3f32..10f32), rng.gen_range(3f32..10f32)];
        tree.insert(&point);
        points.push(point);
    }

    while !points.is_empty() {
        let remove_index = rng.gen_range(0..points.len());
        tree.delete(&points.remove(remove_index));
        getchar::getchar();

        println!("{:?}", tree);
        let mut canvas = Canvas::new(1000, 1000);
        draw_node(&tree.root, &mut canvas, &mut LevelColor::default());
        println!("drawn");
        std::fs::copy("plot.svg", "plot_old.svg");
        render::save(&canvas, "plot.svg", SvgRenderer::new()).unwrap();
    }
}

const SCALE: f32 = 50.0;

struct LevelColor {
    level: usize,
    colors: Vec<RGB>,
}

impl LevelColor {
    pub fn new() -> Self {
        Self {
            level: 0,
            colors: vec![
                RGB { r: 0, g: 0, b: 0 },
                RGB { r: 255, g: 0, b: 0 },
                RGB { r: 0, g: 255, b: 0 },
                RGB { r: 0, g: 0, b: 255 },
                RGB {
                    r: 255,
                    g: 255,
                    b: 0,
                },
                RGB {
                    r: 0,
                    g: 255,
                    b: 255,
                },
                RGB {
                    r: 255,
                    g: 0,
                    b: 255,
                },
            ],
        }
    }

    pub fn get(&mut self) -> RGB {
        if self.colors.len() <= self.level {
            self.colors.push(Color::random())
        }
        self.colors[self.level].clone()
    }
    pub fn inc(&mut self) {
        self.level += 1;
    }
    pub fn dec(&mut self) {
        assert!(self.level > 0);
        self.level -= 1;
    }
}

impl Default for LevelColor {
    fn default() -> Self {
        Self::new()
    }
}

fn draw_node<const K: usize, const M: usize>(
    node: &SsNode<K, M>,
    canvas: &mut Canvas,
    level_color: &mut LevelColor,
) {
    let color = Color::random();
    let centroid_x = node.centroid[0] * SCALE;
    let centroid_y = node.centroid[1] * SCALE;
    let circle = Drawing::new()
        .with_shape(Shape::Circle {
            radius: (node.radius * SCALE) as u32,
        })
        .with_xy(centroid_x, centroid_y)
        .with_style(Style::stroked(1, color));

    canvas.display_list.add(circle);

    let center_cross = Drawing::new()
        .with_shape(
            LineBuilder::new(centroid_x, centroid_y)
                .line_to(centroid_x - 5.0, centroid_y - 5.0)
                .line_to(centroid_x + 5.0, centroid_y - 5.0)
                .line_to(centroid_x - 5.0, centroid_y + 5.0)
                .line_to(centroid_x + 5.0, centroid_y + 5.0)
                .line_to(centroid_x, centroid_y)
                .build(),
        )
        .with_style(Style::stroked(1, color));
    canvas.display_list.add(center_cross);
    match &node.links {
        aadsia::sstree::SsNodeLinks::Inner(nodes) => {
            level_color.inc();
            for node in nodes.iter() {
                draw_node(node, canvas, level_color);
            }
            level_color.dec();
        }
        aadsia::sstree::SsNodeLinks::Leaf(points) => {
            if !false {
                for point in points.iter() {
                    let point_drawing = Drawing::new()
                        .with_shape(Shape::Circle { radius: 2 })
                        .with_xy(point[0] * SCALE, point[1] * SCALE)
                        .with_style(Style::stroked(2, color));
                    canvas.display_list.add(point_drawing);
                }
            }
        }
    }
}
