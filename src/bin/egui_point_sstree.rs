use std::time::Instant;

use aadsia::point_sstree::SsTree;
use eframe::epi;
use egui::{emath, util::cache::CacheTrait, Color32, Frame, Pos2, Rect, Sense, Shape, Stroke};
use rand::Rng;

struct Select {
    center: Pos2,
    radius: f32,
}

fn pos2_to_array(pos: &Pos2) -> [f32; 2] {
    [pos.x, pos.y]
}

impl Select {
    pub fn new(center: Pos2) -> Self {
        Self {
            center,
            radius: 0.0,
        }
    }

    pub fn update(&mut self, pos: Pos2) {
        let d = pos - self.center;
        self.radius = d.length();
    }
}

const M: usize = 32;
const LOWER_M: usize = 16;

#[derive(Default)]
struct MyEguiApp {
    stroke: Stroke,
    shapes: Vec<Shape>,

    tree: SsTree<2, M>,
    max_depth: usize,
    draw_points: bool,
    select: bool,
    select_tool: Option<Select>,
    delete: bool,
}

impl epi::App for MyEguiApp {
    fn update(&mut self, ctx: &egui::CtxRef, frame: &mut epi::Frame<'_>) {
        egui::CentralPanel::default().show(ctx, |ui| {
            let res = ui.add(egui::Slider::new(&mut self.max_depth, 0..=6));
            let res2 = ui.add(egui::Checkbox::new(&mut self.draw_points, "points"));
            let res_select = ui.add(egui::Checkbox::new(&mut self.select, "select"));
            let res_delete = ui.add(egui::Checkbox::new(&mut self.delete, "delete"));

            // println!("res: {:?}", res);
            Frame::dark_canvas(ui.style()).show(ui, |ui| {
                self.ui_canvas(ui, res.changed() || res2.changed());
            });
        });
    }

    fn name(&self) -> &str {
        "SS-Tree"
    }
}

impl MyEguiApp {
    pub fn ui_canvas(&mut self, ui: &mut egui::Ui, mut changed: bool) -> egui::Response {
        // self.ui_content(ui);
        ui.heading("SS-Tree");

        let (mut response, painter) =
            ui.allocate_painter(ui.available_size_before_wrap(), Sense::drag());
        let to_screen = emath::RectTransform::from_to(
            Rect::from_min_size(Pos2::ZERO, response.rect.square_proportions()),
            response.rect,
        );

        if self.select {
            if response.drag_started() {
                self.select_tool = Some(Select::new(
                    response
                        .interact_pointer_pos()
                        .expect("missing pointer pos in drag start"),
                ));
            } else if response.drag_released() {
                self.select_tool = None;
            } else if response.dragged() {
                if let Some(select_tool) = self.select_tool.as_mut() {
                    select_tool.update(
                        response
                            .interact_pointer_pos()
                            .expect("missing pointer pos in drag"),
                    );
                }
            }
        } else if self.delete {
            if response.dragged() {
                let start = Instant::now();
                let mut selected = Vec::new();
                self.tree.points_within_radius(
                    &pos2_to_array(&response.interact_pointer_pos().unwrap()),
                    100.0,
                    &mut selected,
                );
                println!("selected: {} in {:?}", selected.len(), start.elapsed());
                let start = Instant::now();
                changed = !selected.is_empty();
                for point in selected {
                    self.tree.delete(&point);
                }
                println!("deleted: {:?}", start.elapsed());
            }
        } else {
            let from_screen = to_screen.inverse();
            if let Some(pointer_pos) = response.interact_pointer_pos() {
                let canvas_pos = from_screen * pointer_pos;
                println!("interact: {:?}", canvas_pos);
                self.tree.insert(&[pointer_pos.x, pointer_pos.y]);
                changed = true;
                // self.shapes
                //     .push(egui::Shape::circle_stroke(pointer_pos, 10.0, self.stroke));
                // if current_line.last() != Some(&canvas_pos) {
                //     current_line.push(canvas_pos);
                //     response.mark_changed();
                // }
                // response.mark_changed();
            }
        }

        if changed {
            self.shapes.clear();
            draw_tree(
                &mut self.shapes,
                &self.tree.root,
                0,
                self.max_depth,
                self.draw_points,
            );
        }

        painter.extend(self.shapes.clone());

        if let Some(select_tool) = self.select_tool.as_ref() {
            painter.add(egui::Shape::circle_stroke(
                select_tool.center,
                select_tool.radius,
                Stroke::new(1.0, Color32::WHITE),
            ));

            let mut selected = Vec::new();

            let start = Instant::now();
            self.tree.points_within_radius(
                &[select_tool.center.x, select_tool.center.y],
                select_tool.radius,
                &mut selected,
            );
            println!("selected: {} {:?}", selected.len(), start.elapsed());

            painter.extend(
                selected
                    .iter()
                    .map(|p| {
                        let center = Pos2::new(p[0], p[1]);
                        egui::Shape::circle_filled(center, 1.5, Color32::BLUE)
                    })
                    .collect(),
            );
        }

        response
    }
}

const COLORS: [Color32; 9] = [
    Color32::RED,
    Color32::GREEN,
    Color32::BLUE,
    Color32::LIGHT_RED,
    Color32::LIGHT_GREEN,
    Color32::LIGHT_BLUE,
    Color32::DARK_RED,
    Color32::DARK_GREEN,
    Color32::DARK_BLUE,
];

fn draw_tree<const K: usize, const M: usize>(
    shapes: &mut Vec<Shape>,
    node: &aadsia::point_sstree::SsNode<K, M>,
    level: usize,
    max_level: usize,
    draw_points: bool,
) {
    if level > max_level {
        return;
    }
    let circle = egui::Shape::circle_stroke(
        Pos2::new(node.centroid[0], node.centroid[1]),
        node.radius,
        Stroke::new(1.0, COLORS[level]),
    );
    shapes.push(circle);
    // let circle = Drawing::new()
    //     .with_shape(Shape::Circle {
    //         radius: (node.radius * SCALE) as u32,
    //     })
    //     .with_xy(node.centroid[0] * SCALE, node.centroid[1] * SCALE)
    //     .with_style(Style::stroked(1, level_color.get()));

    // canvas.display_list.add(circle);

    match &node.links {
        aadsia::point_sstree::SsNodeLinks::Inner(nodes) => {
            // level_color.inc();
            for node in nodes.iter() {
                draw_tree(shapes, node, level + 1, max_level, draw_points);
            }
            // level_color.dec();
        }
        aadsia::point_sstree::SsNodeLinks::Leaf(points) => {
            if draw_points {
                for point in points.iter() {
                    let point = egui::Shape::circle_filled(
                        Pos2::new(point[0], point[1]),
                        1.0,
                        Color32::WHITE,
                    );
                    shapes.push(point);
                    // let point_drawing = Drawing::new()
                    //     .with_shape(Shape::Circle { radius: 1 })
                    //     .with_xy(point[0] * SCALE, point[1] * SCALE)
                    //     .with_style(Style::stroked(1, Color::black()));
                    // canvas.display_list.add(point_drawing);
                }
            }
        }
    }
}

fn main() {
    let mut tree = SsTree::new(LOWER_M);
    let mut rng = rand::thread_rng();

    // for _ in 0..1000000 {
    //     tree.insert(&[rng.gen_range(200.0..600.0), rng.gen_range(200.0..600.0)]);
    // }

    let app = MyEguiApp {
        stroke: Stroke::new(1.0, Color32::LIGHT_BLUE),
        shapes: Vec::new(),
        tree: tree,
        max_depth: 2,
        draw_points: true,
        select: false,
        select_tool: None,
        delete: false,
    };
    let native_options = eframe::NativeOptions::default();
    eframe::run_native(Box::new(app), native_options);
}
