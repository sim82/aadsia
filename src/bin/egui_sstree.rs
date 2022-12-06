use std::time::Instant;

use aadsia::sstree::{Dimindex, Distance, Element, SsTree};
use eframe::epi;
use egui::{emath, Color32, Frame, Pos2, Rect, Sense, Shape, Stroke};
use rand::Rng;

struct Select {
    center: Pos2,
    radius: f32,
}

#[derive(Debug)]
struct Drag {
    path: Vec<u8>,
    last_pos: Pos2,
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

impl Drag {
    fn update<P, K: Distance + Dimindex + PartialEq + Default, const M: usize>(
        &mut self,
        pos: Pos2,
        tree: &mut SsTree<P, K, M>,
    ) {
        let mut element = tree.remove_by_path(&self.path);

        let d = pos - self.last_pos;
        element.center[0] += d.x;
        element.center[1] += d.y;
        self.last_pos = pos;
        self.path = tree.insert_get_path(element);
    }
}

const M: usize = 32;
const LOWER_M: usize = 16;

#[derive(Default, PartialEq)]
enum Mode {
    #[default]
    Draw,
    Select,
    Delete,
    Drag,
}

#[derive(Default)]
struct MyEguiApp {
    shapes: Vec<Shape>,

    tree: SsTree<u32, [f32; 2], M>,

    mode: Mode,
    max_depth: usize,
    draw_points: bool,
    select_tool: Option<Select>,
    smear: bool,
    delete_radius: f32,
    insert_radius: f32,
    insert_count: u32,
    drag_tool: Option<Drag>,
}

impl epi::App for MyEguiApp {
    fn update(&mut self, ctx: &egui::CtxRef, _frame: &mut epi::Frame<'_>) {
        egui::CentralPanel::default().show(ctx, |ui| {
            let res = ui.add(egui::Slider::new(&mut self.max_depth, 0..=10));
            let res2 = ui.add(egui::Checkbox::new(&mut self.draw_points, "points"));

            ui.radio_value(&mut self.mode, Mode::Draw, "draw");
            ui.radio_value(&mut self.mode, Mode::Select, "select");
            ui.radio_value(&mut self.mode, Mode::Delete, "delete");
            ui.radio_value(&mut self.mode, Mode::Drag, "drag");

            // ui.add(egui::Checkbox::new(&mut self.select, "select"));
            // ui.add(egui::Checkbox::new(&mut self.delete, "delete"));
            // ui.add(egui::Checkbox::new(&mut self.drag, "drag"));
            ui.add(egui::Checkbox::new(&mut self.smear, "smear"));
            ui.add(egui::Slider::new(&mut self.insert_radius, 1.0..=20.0));
            ui.add(egui::Slider::new(&mut self.delete_radius, 5.0..=100.0).text("delete radius"));

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

        let (response, painter) =
            ui.allocate_painter(ui.available_size_before_wrap(), Sense::drag());
        let to_screen = emath::RectTransform::from_to(
            Rect::from_min_size(Pos2::ZERO, response.rect.square_proportions()),
            response.rect,
        );

        match self.mode {
            Mode::Select => {
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
            }
            Mode::Delete if response.dragged() => {
                let start = Instant::now();

                let mut selected = Vec::new();

                self.tree.points_within_radius(
                    &pos2_to_array(&response.interact_pointer_pos().unwrap()),
                    self.delete_radius,
                    &mut selected,
                );
                painter.add(egui::Shape::circle_stroke(
                    response.interact_pointer_pos().unwrap(),
                    self.delete_radius,
                    Stroke::new(1.0, Color32::WHITE),
                ));

                println!("selected: {} in {:?}", selected.len(), start.elapsed());
                let start = Instant::now();
                changed = !selected.is_empty();
                let centers = selected.drain(..).map(|e| e.center).collect::<Vec<_>>();
                for point in centers {
                    self.tree.delete(&point); // FIXME: we probably want to delete by identity
                }
                println!("deleted: {:?}", start.elapsed());
            }
            Mode::Drag => {
                if response.drag_started() {
                    let mut selected = Vec::new();

                    let pos2 = response.interact_pointer_pos().unwrap();
                    self.tree
                        .paths_within_radius(&pos2_to_array(&pos2), 1.0, &mut selected);
                    println!("drag start {:?}", selected);

                    if let Some(path) = selected.pop() {
                        self.drag_tool = Some(Drag {
                            path,
                            last_pos: pos2,
                        });
                        // self.tree.remove_by_path(&path);
                        changed = true;
                    }
                } else if response.drag_released() {
                    self.drag_tool = None;
                    changed = true;
                } else if response.dragged() {
                    if let Some(drag_tool) = self.drag_tool.as_mut() {
                        drag_tool.update(
                            response
                                .interact_pointer_pos()
                                .expect("missing pointer pos in drag"),
                            &mut self.tree,
                        );
                    }
                    changed = true;
                    println!("drag update {:?}", self.drag_tool);
                }
            }
            Mode::Draw => {
                if response.clicked() || self.smear {
                    let from_screen = to_screen.inverse();

                    if let Some(pointer_pos) = response.interact_pointer_pos() {
                        let canvas_pos = from_screen * pointer_pos;
                        println!("interact: {:?}", canvas_pos);
                        let path = self.tree.insert_get_path(Element::new(
                            [pointer_pos.x, pointer_pos.y],
                            self.insert_radius,
                            self.insert_count,
                        ));
                        self.insert_count += 1;
                        changed = true;
                        println!("insert at: {:?}", path);
                        // self.shapes
                        //     .push(egui::Shape::circle_stroke(pointer_pos, 10.0, self.stroke));
                        // if current_line.last() != Some(&canvas_pos) {
                        //     current_line.push(canvas_pos);
                        //     response.mark_changed();
                        // }
                        // response.mark_changed();
                    }
                }
            }
            _ => (),
        }
        if changed {
            self.shapes.clear();
            draw_tree(
                &mut self.shapes,
                &self.tree.root,
                self.max_depth,
                self.draw_points,
            );
        }

        painter.extend(self.shapes.clone());
        // for i in 0..100 {
        //     tree.insert(Element::new(
        //         [rng.gen_range(200.0..600.0), rng.gen_range(200.0..600.0)],
        //         2.0,
        //         i,
        //     ));
        // }

        if let Some(select_tool) = self.select_tool.as_ref() {
            painter.add(egui::Shape::circle_stroke(
                select_tool.center,
                select_tool.radius,
                Stroke::new(1.0, Color32::WHITE),
            ));

            let mut selected = Vec::new();

            let start = Instant::now();
            if !false {
                self.tree.points_within_radius(
                    &[select_tool.center.x, select_tool.center.y],
                    select_tool.radius,
                    &mut selected,
                );
            } else {
                let mut paths = Vec::new();
                self.tree.paths_within_radius(
                    &[select_tool.center.x, select_tool.center.y],
                    select_tool.radius,
                    &mut paths,
                );

                for path in paths {
                    selected.push(self.tree.get_by_path(&path))
                }
            }
            println!("selected: {} {:?}", selected.len(), start.elapsed());

            painter.extend(
                selected
                    .iter()
                    .map(|p| {
                        let center = Pos2::new(p.center[0], p.center[1]);
                        egui::Shape::circle_filled(center, p.radius + 0.5, Color32::BLUE)
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

fn draw_tree<const M: usize>(
    shapes: &mut Vec<Shape>,
    node: &aadsia::sstree::SsNode<u32, [f32; 2], M>,
    max_level: usize,
    draw_points: bool,
) {
    let mut stack = Vec::new();
    stack.push((node, 0));
    let mut element_color = 0;
    while let Some((node, level)) = stack.pop() {
        if level > max_level {
            continue;
        }

        let circle = egui::Shape::circle_stroke(
            Pos2::new(node.centroid[0], node.centroid[1]),
            node.radius,
            Stroke::new(1.0, COLORS[level]),
        );

        shapes.push(circle);
        shapes.push(egui::Shape::circle_stroke(
            Pos2::new(node.centroid[0], node.centroid[1]),
            1.0,
            Stroke::new(1.0, COLORS[level]),
        ));

        match &node.links {
            aadsia::sstree::SsNodeLinks::Inner(nodes) => {
                // level_color.inc();
                for node in nodes.iter() {
                    stack.push((node, level + 1))
                }
                // level_color.dec();
            }
            aadsia::sstree::SsNodeLinks::Leaf(points) => {
                if draw_points {
                    for point in points.iter() {
                        let point = egui::Shape::circle_filled(
                            Pos2::new(point.center[0], point.center[1]),
                            point.radius,
                            COLORS[element_color % COLORS.len()],
                        );
                        shapes.push(point);
                    }
                    element_color += 1;
                }
            }
        }
    }
}

fn main() {
    let mut tree = SsTree::new(LOWER_M);
    let mut rng = rand::thread_rng();

    for i in 0..1000000 {
        tree.insert(Element::new(
            [rng.gen_range(200.0..60000.0), rng.gen_range(200.0..60000.0)],
            5.0,
            i,
        ));
    }

    let app = MyEguiApp {
        shapes: Vec::new(),
        tree,
        mode: Mode::Draw,
        max_depth: 2,
        draw_points: true,
        select_tool: None,
        smear: false,
        insert_radius: 5.0,
        delete_radius: 20.0,
        insert_count: 0,
        drag_tool: None,
    };
    let native_options = eframe::NativeOptions::default();
    eframe::run_native(Box::new(app), native_options);
}
