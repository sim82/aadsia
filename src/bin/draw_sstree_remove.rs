use aadsia::sstree::{SsNode, SsTree};
use draw::{render, Canvas, Color, Drawing, Shape, Style, SvgRenderer, RGB};
use rand::{random, Rng};

fn main() {
    let mut tree = SsTree::<2, 8>::new(4);

    tree.insert(&[1.0, 1.0]);
    tree.insert(&[0.0, 0.0]);
    tree.insert(&[2.0, 2.0]);
    tree.insert(&[3.0, 3.0]);

    println!("{:?}", tree);

    tree.insert(&[4.0, 4.0]);
    tree.insert(&[5.0, 5.0]);
    tree.insert(&[6.0, 6.0]);
    tree.insert(&[7.0, 7.0]);

    println!("{:?}", tree);
    tree.insert(&[9.0, 9.0]);

    tree.delete(&[5.0, 5.0]);
    tree.delete(&[2.0, 2.0]);
    // tree.insert(&[4.0, 4.0]);

    println!("{:?}", tree);
    let mut canvas = Canvas::new(1000, 1000);
    draw_node(&tree.root, &mut canvas, &mut LevelColor::default());
    println!("drawed");
    render::save(&canvas, "plot.svg", SvgRenderer::new()).unwrap();
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
    let circle = Drawing::new()
        .with_shape(Shape::Circle {
            radius: (node.radius * SCALE) as u32,
        })
        .with_xy(node.centroid[0] * SCALE, node.centroid[1] * SCALE)
        .with_style(Style::stroked(1, level_color.get()));

    canvas.display_list.add(circle);

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
                        .with_shape(Shape::Circle { radius: 1 })
                        .with_xy(point[0] * SCALE, point[1] * SCALE)
                        .with_style(Style::stroked(1, Color::black()));
                    canvas.display_list.add(point_drawing);
                }
            }
        }
    }
}
