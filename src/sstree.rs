use arrayvec::ArrayVec;

pub trait Distance {
    fn distance(&self, other: &Self) -> f32;
}

pub trait Zero {
    fn zero() -> Self;
}

pub trait MaxVarianceDirection: Sized {
    fn mean_along_direction(elements: &[&Self], direction_index: usize) -> f32;
    fn variance_along_direction(elements: &[&Self], direction_index: usize) -> f32;
    fn direction_of_max_variance(elements: &[&Self]) -> usize;
    fn centroid(elements: &[&Self]) -> Self;
    fn dim_component(&self, i: usize) -> f32;
}

impl<const K: usize> MaxVarianceDirection for [f32; K] {
    fn dim_component(&self, i: usize) -> f32 {
        self[i]
    }

    fn mean_along_direction(elements: &[&[f32; K]], direction_index: usize) -> f32 {
        assert!(!elements.is_empty());
        let count = elements.len() as f32;
        let sum = elements
            .iter()
            .map(|point| point[direction_index])
            .sum::<f32>();
        sum / count
    }

    fn variance_along_direction(elements: &[&[f32; K]], direction_index: usize) -> f32 {
        assert!(!elements.is_empty());
        let mean = Self::mean_along_direction(elements, direction_index);
        let count = elements.len() as f32;
        elements
            .iter()
            .map(|point| {
                let diff = mean - point[direction_index];

                diff * diff
            })
            .sum::<f32>()
            / count
    }

    fn direction_of_max_variance(elements: &[&[f32; K]]) -> usize {
        let mut max_variance = 0.0;
        let mut direction_index = 0;
        for i in 0..K {
            let variance = Self::variance_along_direction(elements, i);
            if variance > max_variance {
                max_variance = variance;
                direction_index = i;
            }
        }
        direction_index
    }

    fn centroid(elements: &[&[f32; K]]) -> Self {
        let mut centroid = [0f32; K];
        for (i, c) in centroid.iter_mut().enumerate().take(K) {
            *c = Self::mean_along_direction(elements, i);
        }
        centroid
    }
}

#[derive(Debug)]
pub struct Element<P, K> {
    pub center: K,
    pub radius: f32,
    pub payload: P,
}

impl<P, K: PartialEq> PartialEq for Element<P, K> {
    fn eq(&self, other: &Self) -> bool {
        self.center == other.center && self.radius == other.radius
    }
}

impl<P, K: Distance> Element<P, K> {
    pub fn new(center: K, radius: f32, payload: P) -> Self {
        Self {
            center,
            radius,
            payload,
        }
    }
    pub fn intersects_point(&self, target: &K) -> bool {
        self.center.distance(target) <= self.radius
    }
}

#[derive(Debug)]
pub enum SsNodeLinks<P, K: Distance + MaxVarianceDirection + PartialEq, const M: usize> {
    Inner(Box<ArrayVec<SsNode<P, K, M>, M>>),
    Leaf(Box<ArrayVec<Element<P, K>, M>>),
}

#[derive(Debug)]
pub struct SsNode<P, K: Distance + MaxVarianceDirection + PartialEq, const M: usize> {
    pub centroid: K,
    pub radius: f32,
    pub links: SsNodeLinks<P, K, M>,
}

impl<const K: usize> Distance for [f32; K] {
    fn distance(&self, p2: &[f32; K]) -> f32 {
        self.iter()
            .zip(p2.iter())
            .map(|(c1, c2)| (c1 - c2) * (c1 - c2))
            .sum::<f32>()
            .sqrt()
    }
}
impl<const K: usize> Zero for [f32; K] {
    fn zero() -> Self {
        [0f32; K]
    }
}

#[test]
fn test_distance() {
    assert_eq!([0.0, 0.0].distance(&[1.0, 1.0]), 2.0f32.sqrt());
    assert_eq!([-10.0, 1.0].distance(&[10.0, 1.0]), 20.0f32);
    assert_eq!([1000.0, -1000.0].distance(&[1000.0, 2000.0]), 3000.0);
}

impl<P, K: Zero + MaxVarianceDirection + Distance + PartialEq, const M: usize> SsNode<P, K, M> {
    pub fn from_elements(elements: ArrayVec<Element<P, K>, M>) -> Self {
        let (centroid, radius) = leaf::centroid_and_radius::<P, K, M>(&elements);
        Self {
            centroid,
            radius,
            links: SsNodeLinks::Leaf(Box::new(elements)),
        }
    }

    pub fn from_nodes(nodes: ArrayVec<Self, M>) -> Self {
        let (centroid, radius) = inner::centroid_and_radius(&nodes);
        Self {
            centroid,
            radius,
            links: SsNodeLinks::Inner(Box::new(nodes)),
        }
    }

    pub fn intersects_point(&self, target: &K) -> bool {
        self.centroid.distance(target) <= self.radius
    }

    pub fn search(&self, target: &K) -> Option<&Self> {
        match &self.links {
            SsNodeLinks::Inner(children) => {
                children.iter().find(|node| node.intersects_point(target))
            }
            SsNodeLinks::Leaf(points) => {
                if points.iter().any(|x| x.intersects_point(target)) {
                    Some(self)
                } else {
                    None
                }
            }
        }
    }

    pub fn search_parent_leaf(&self, target: &K) -> &Self {
        match &self.links {
            SsNodeLinks::Inner(children) => {
                let child = find_closest_child(children, target);
                child.search_parent_leaf(target)
            }
            SsNodeLinks::Leaf(_) => self,
        }
    }

    pub fn update_bounding_envelope(&mut self) {
        let (centroid, radius) = match &self.links {
            SsNodeLinks::Inner(nodes) => inner::centroid_and_radius(nodes),
            SsNodeLinks::Leaf(points) => leaf::centroid_and_radius::<P, K, M>(points),
        };
        self.centroid = centroid;
        self.radius = radius;
    }
    pub fn insert(&mut self, element: Element<P, K>, m: usize) -> Option<(Self, Self)> {
        match &mut self.links {
            SsNodeLinks::Leaf(points) => {
                if points.iter().any(|p| *p == element) {
                    return None;
                }

                if points.len() < M {
                    points.push(element);
                    self.update_bounding_envelope();
                    return None;
                } else {
                    let mut nodes_to_split = points
                        .drain(..)
                        .chain(std::iter::once(element))
                        .collect::<Vec<_>>();

                    let split_index = leaf::find_split_index::<P, K, M>(&mut nodes_to_split, m);
                    let points2: ArrayVec<_, M> = nodes_to_split.drain(split_index..).collect();
                    let (centroid2, radius2) = leaf::centroid_and_radius::<P, K, M>(&points2);

                    let points1: ArrayVec<_, M> = nodes_to_split.drain(..split_index).collect();
                    let (centroid1, radius1) = leaf::centroid_and_radius::<P, K, M>(&points1);

                    let new_node1 = Self {
                        centroid: centroid1,
                        radius: radius1,
                        links: SsNodeLinks::Leaf(Box::new(points1)),
                    };
                    let new_node2 = Self {
                        centroid: centroid2,
                        radius: radius2,
                        links: SsNodeLinks::Leaf(Box::new(points2)),
                    };
                    return Some((new_node1, new_node2));
                }
            }

            SsNodeLinks::Inner(children) => {
                // TODO: check if using element.centroid also works for non-point elements
                let closest_child_index = find_closest_child_index(children, &element.center);
                if let Some((new_child_1, new_child_2)) =
                    children[closest_child_index].insert(element, m)
                {
                    children.remove(closest_child_index);

                    if children.len() < M - 1 {
                        children.push(new_child_1);
                        children.push(new_child_2);
                    } else {
                        let mut nodes_to_split: Vec<_> = children
                            .drain(..)
                            .chain(std::iter::once(new_child_1))
                            .chain(std::iter::once(new_child_2))
                            .collect();

                        let split_index = inner::find_split_index(&mut nodes_to_split, m);

                        let points2: ArrayVec<_, M> = nodes_to_split.drain(split_index..).collect();
                        let (centroid2, radius2) = inner::centroid_and_radius(&points2);

                        let points1: ArrayVec<_, M> = nodes_to_split.drain(..split_index).collect();
                        let (centroid1, radius1) = inner::centroid_and_radius(&points1);

                        let new_node1 = Self {
                            centroid: centroid1,
                            radius: radius1,
                            links: SsNodeLinks::Inner(Box::new(points1)),
                        };
                        let new_node2 = Self {
                            centroid: centroid2,
                            radius: radius2,
                            links: SsNodeLinks::Inner(Box::new(points2)),
                        };
                        return Some((new_node1, new_node2));
                    }
                } else {
                    self.update_bounding_envelope();
                }
            }
        }
        None
    }

    pub fn delete(&mut self, target: &K, m: usize) -> (bool, bool) {
        match &mut self.links {
            SsNodeLinks::Leaf(elements) => {
                if let Some((i, _)) = elements
                    .iter()
                    .enumerate()
                    .find(|(_, p)| p.intersects_point(target))
                {
                    elements.remove(i);
                    let num_elemens = elements.len();
                    self.update_bounding_envelope();
                    (true, num_elemens < m)
                } else {
                    (false, false)
                }
            }
            SsNodeLinks::Inner(nodes) => {
                let mut node_to_fix_index = None;
                let mut deleted = false;
                for (i, child_node) in nodes.iter_mut().enumerate() {
                    if child_node.intersects_point(target) {
                        let res = child_node.delete(target, m);
                        deleted = res.0;
                        let violates_invariants = res.1;
                        // println!("{:?} {:?}", deleted, violates_invariants);
                        if violates_invariants {
                            node_to_fix_index = Some(i);
                        }
                        if deleted {
                            break;
                        }
                    }
                }
                match node_to_fix_index {
                    None => {
                        if deleted {
                            self.update_bounding_envelope();
                        }
                        (deleted, false)
                    }

                    Some(node_to_fix) => {
                        if let Some(sibling_to_borrow_from) =
                            inner::find_sibling_to_borrow_from(nodes, node_to_fix, m)
                        {
                            inner::borrow_from_sibling(nodes, node_to_fix, sibling_to_borrow_from);
                        } else if let Some(sibling_to_merge_to) =
                            inner::find_sibling_to_merge_to(nodes, node_to_fix, m)
                        {
                            // no sibling to borrow from -> merge
                            inner::merge_siblings(nodes, node_to_fix, sibling_to_merge_to);
                        }
                        //  else {
                        //     panic!("cannot borrow / merge");
                        // }
                        (true, nodes.len() < m)
                    }
                }
            }
        }
    }

    pub fn count_nodes(&self) -> (usize, usize) {
        match &self.links {
            SsNodeLinks::Inner(nodes) => nodes.iter().fold((0, 1), |(a_points, a_nodes), n| {
                let (points, nodes) = n.count_nodes();
                (a_points + points, a_nodes + nodes)
            }),
            SsNodeLinks::Leaf(points) => (points.len(), 1),
        }
    }
    pub fn elements_within_radius<'a>(
        &'a self,
        center: &K,
        radius: f32,
        out: &mut Vec<&'a Element<P, K>>,
    ) {
        match &self.links {
            SsNodeLinks::Leaf(points) => {
                for point in points.iter() {
                    if point.center.distance(center) < (radius + point.radius) {
                        out.push(point);
                    }
                }
            }
            SsNodeLinks::Inner(nodes) => {
                for child in nodes.iter() {
                    if child.centroid.distance(center) <= radius + child.radius {
                        child.elements_within_radius(center, radius, out);
                    }
                }
            }
        }
    }
    // function pointsWithinRegion(node, region)
    //   points ← []
    //   if node.leaf then
    //     for point in node.points do
    //       if region.intersectsPoint(point) then
    //         points.insert(point)
    //   else
    //    for child in node.children do
    //       if region.intersectsNode(child) then
    //         points.insertAll(pointsWithinRegion(child, region))
    //   return points
}

fn find_closest_child<'a, P, K: Distance + MaxVarianceDirection + PartialEq, const M: usize>(
    children: &'a [SsNode<P, K, M>],
    target: &K,
) -> &'a SsNode<P, K, M> {
    // children
    //     .iter()
    //     .min_by_key(|a| distance(&a.centroid, target))
    //     .as_ref()
    //     .unwrap()

    let mut min_dist = f32::MAX;
    let mut cur_min = None;
    for child in children {
        let d = child.centroid.distance(target);
        if d < min_dist {
            min_dist = d;
            cur_min = Some(child);
        }
    }
    cur_min.unwrap()
}
fn find_closest_child_index<P, K: Distance + MaxVarianceDirection + PartialEq, const M: usize>(
    children: &[SsNode<P, K, M>],
    target: &K,
) -> usize {
    let mut min_dist = f32::MAX;
    let mut cur_min = None;
    for (i, child) in children.iter().enumerate() {
        let d = child.centroid.distance(target);
        if d < min_dist {
            min_dist = d;
            cur_min = Some(i);
        }
    }
    cur_min.unwrap()
}

#[derive(Debug)]
pub struct SsTree<P, K: Distance + MaxVarianceDirection + PartialEq, const M: usize> {
    pub root: SsNode<P, K, M>,
    height: usize,
    m: usize,
}

impl<P, K: Zero + Distance + MaxVarianceDirection + PartialEq, const M: usize> SsTree<P, K, M> {
    pub fn new(m: usize) -> Self {
        Self {
            root: SsNode {
                centroid: K::zero(),
                radius: 0f32,
                links: SsNodeLinks::Leaf(Box::new(ArrayVec::new())),
            },
            height: 1,
            m,
        }
    }

    pub fn insert(&mut self, element: Element<P, K>) {
        if let Some((new_child_1, new_child_2)) = self.root.insert(element, self.m) {
            let mut nodes = ArrayVec::<_, M>::new();
            nodes.push(new_child_1);
            nodes.push(new_child_2);
            let (centroid, radius) = inner::centroid_and_radius(&nodes);
            self.root = SsNode {
                centroid,
                radius,
                links: SsNodeLinks::Inner(Box::new(nodes)),
            };
            self.height += 1;
        }
    }
    #[allow(clippy::overly_complex_bool_expr)]
    pub fn delete(&mut self, point: &K) {
        let (_deleted, violiates_invariant) = self.root.delete(point, self.m);

        if false && violiates_invariant {
            println!("violate invariant at root!");
            match &mut self.root.links {
                SsNodeLinks::Leaf(leafs) => println!("leafs: {}", leafs.len()),
                SsNodeLinks::Inner(nodes) => {
                    println!("nodes: {}", nodes.len());
                    self.root = nodes.pop().unwrap();
                }
            }
        }
    }

    pub fn get_height(&self) -> usize {
        self.height
    }
    pub fn get_fill_factor(&self) -> f32 {
        let (num_points, num_nodes) = self.root.count_nodes();
        num_points as f32 / num_nodes as f32
    }

    pub fn points_within_radius<'a>(
        &'a self,
        center: &K,
        radius: f32,
        out: &mut Vec<&'a Element<P, K>>,
    ) {
        self.root.elements_within_radius(center, radius, out);
    }
}

impl<P, K: Distance + MaxVarianceDirection + Zero + PartialEq, const M: usize> Default
    for SsTree<P, K, M>
{
    fn default() -> Self {
        Self::new(M / 2)
    }
}

mod leaf {
    use arrayvec::ArrayVec;

    use super::{Distance, Element, MaxVarianceDirection, Zero};
    // pub fn mean_along_direction<P, K>(elements: &[Element<P, K>], direction_index: usize) -> f32 {
    //     assert!(!elements.is_empty());
    //     let count = elements.len() as f32;
    //     let sum = elements
    //         .iter()
    //         .map(|point| point.center[direction_index])
    //         .sum::<f32>();
    //     sum / count
    // }

    // pub fn variance_along_direction<P, K>(
    //     elements: &[Element<P, K>],
    //     direction_index: usize,
    // ) -> f32 {
    //     assert!(!elements.is_empty());
    //     let mean = mean_along_direction(elements, direction_index);
    //     let count = elements.len() as f32;
    //     elements
    //         .iter()
    //         .map(|point| {
    //             let diff = mean - point.center[direction_index];

    //             diff * diff
    //         })
    //         .sum::<f32>()
    //         / count
    // }

    // pub fn direction_of_max_variance<P, K>(elements: &[Element<P, K>]) -> usize {
    //     let mut max_variance = 0.0;
    //     let mut direction_index = 0;
    //     for i in 0..K {
    //         let variance = variance_along_direction(elements, i);
    //         if variance > max_variance {
    //             max_variance = variance;
    //             direction_index = i;
    //         }
    //     }
    //     direction_index
    // }

    pub fn find_split_index<P, K: MaxVarianceDirection + Distance, const M: usize>(
        elements: &mut [Element<P, K>],
        m: usize,
    ) -> usize {
        let coordinate_index = {
            // FIXME: use ArrayVec one generic const exprs are stable
            let keys = elements.iter().map(|e| &e.center).collect::<Vec<_>>();

            K::direction_of_max_variance(&keys)
        };
        elements.sort_by(|p1, p2| {
            p1.center
                .dim_component(coordinate_index)
                .partial_cmp(&p2.center.dim_component(coordinate_index))
                .unwrap()
        });

        // FIXME: use ArrayVec one generic const exprs are stable
        let keys = elements.iter().map(|e| &e.center).collect::<Vec<_>>();
        let mut min_variance = f32::INFINITY;
        let mut split_index = m;
        for i in m..=(elements.len() - m) {
            let variance1 = K::variance_along_direction(&keys[..i], coordinate_index);
            let variance2 = K::variance_along_direction(&keys[i..], coordinate_index);
            let variance = variance1 + variance2;
            if variance < min_variance {
                min_variance = variance;
                split_index = i;
            }
        }
        split_index
    }

    pub fn centroid_and_radius<P, K: Zero + MaxVarianceDirection + Distance, const M: usize>(
        elements: &[Element<P, K>],
    ) -> (K, f32) {
        let keys = elements
            .iter()
            .map(|e| &e.center)
            .collect::<ArrayVec<_, M>>();

        let centroid = K::centroid(&keys);

        // let centroid = mean_along_all_directions(nodes);

        let radius = elements
            .iter()
            .map(|node| centroid.distance(&node.center) + node.radius)
            .max_by(|d1, d2| d1.partial_cmp(d2).unwrap())
            .unwrap();
        (centroid, radius)
    }
}
mod inner {
    use arrayvec::ArrayVec;

    use super::{Distance, MaxVarianceDirection, SsNode, SsNodeLinks, Zero};

    // pub fn find_split_index<P, const K: usize, const M: usize>(
    //     nodes: &mut [SsNode<P, K, M>],
    //     m: usize,
    // ) -> usize {
    //     let coordinate_index = direction_of_max_variance(nodes);
    //     nodes.sort_by(|p1, p2| {
    //         p1.centroid[coordinate_index]
    //             .partial_cmp(&p2.centroid[coordinate_index])
    //             .unwrap()
    //     });
    //     let mut min_variance = f32::INFINITY;
    //     let mut split_index = m;
    //     for i in m..=(nodes.len() - m) {
    //         let variance1 = variance_along_direction(&nodes[..i], coordinate_index);
    //         let variance2 = variance_along_direction(&nodes[i..], coordinate_index);
    //         let variance = variance1 + variance2;
    //         if variance < min_variance {
    //             min_variance = variance;
    //             split_index = i;
    //         }
    //     }
    //     split_index
    // }

    // pub fn centroid_and_radius<P, const K: usize, const M: usize>(
    //     nodes: &[SsNode<P, K, M>],
    // ) -> ([f32; K], f32) {
    //     let mut centroid = [0f32; K];
    //     for (i, c) in centroid.iter_mut().enumerate().take(K) {
    //         *c = mean_along_direction(nodes, i);
    //     }

    //     // let centroid = mean_along_all_directions(nodes);

    //     let radius = nodes
    //         .iter()
    //         .map(|node| distance(&centroid, &node.centroid) + node.radius)
    //         .max_by(|d1, d2| d1.partial_cmp(d2).unwrap())
    //         .unwrap();
    //     (centroid, radius)
    // }

    // pub fn mean_along_direction<P, const K: usize, const M: usize>(
    //     points: &[SsNode<P, K, M>],
    //     direction_index: usize,
    // ) -> f32 {
    //     assert!(!points.is_empty());
    //     let count = points.len() as f32;
    //     let sum = points
    //         .iter()
    //         .map(|point| point.centroid[direction_index])
    //         .sum::<f32>();
    //     sum / count
    // }

    // // pub fn mean_along_all_directions<const K: usize, const M: usize>(
    // //     points: &[Box<SsNode<K, M>>],
    // // ) -> [f32; K] {
    // //     assert!(!points.is_empty());
    // //     let count = points.len() as f32;
    // //     let mut sum = points.iter().fold([0f32; K], |a, node| {
    // //         let mut a = a.clone();
    // //         a.iter_mut()
    // //             .zip(node.centroid.iter())
    // //             .for_each(|(a, c)| *a += *c);
    // //         a
    // //     });
    // //     sum.iter_mut().for_each(|sum| *sum /= count);
    // //     sum
    // // }

    // pub fn variance_along_direction<P, const K: usize, const M: usize>(
    //     points: &[SsNode<P, K, M>],
    //     direction_index: usize,
    // ) -> f32 {
    //     assert!(!points.is_empty());
    //     let mean = mean_along_direction(points, direction_index);
    //     let count = points.len() as f32;
    //     points
    //         .iter()
    //         .map(|point| {
    //             let diff = mean - point.centroid[direction_index];

    //             diff * diff
    //         })
    //         .sum::<f32>()
    //         / count
    // }

    // pub fn direction_of_max_variance<P, const K: usize, const M: usize>(
    //     nodes: &[SsNode<P, K, M>],
    // ) -> usize {
    //     let mut max_variance = 0.0;
    //     let mut direction_index = 0;
    //     for i in 0..K {
    //         let variance = variance_along_direction(nodes, i);
    //         if variance > max_variance {
    //             max_variance = variance;
    //             direction_index = i;
    //         }
    //     }
    //     direction_index
    // }
    pub fn centroid_and_radius<
        P,
        K: MaxVarianceDirection + Distance + PartialEq,
        const M: usize,
    >(
        nodes: &[SsNode<P, K, M>],
    ) -> (K, f32) {
        let keys = nodes
            .iter()
            .map(|n| &n.centroid)
            .collect::<ArrayVec<_, M>>();

        let centroid = K::centroid(&keys);
        let radius = nodes
            .iter()
            .map(|node| centroid.distance(&node.centroid) + node.radius)
            .max_by(|d1, d2| d1.partial_cmp(d2).unwrap())
            .unwrap();
        (centroid, radius)
    }

    pub fn find_split_index<P, K: Distance + MaxVarianceDirection + PartialEq, const M: usize>(
        nodes: &mut [SsNode<P, K, M>],
        m: usize,
    ) -> usize {
        // FIXME: use ArrayVec when generco const exprs are stable
        let coordinate_index = {
            let keys = nodes.iter().map(|n| &n.centroid).collect::<Vec<_>>();

            K::direction_of_max_variance(&keys)
        };
        nodes.sort_by(|p1, p2| {
            p1.centroid
                .dim_component(coordinate_index)
                .partial_cmp(&p2.centroid.dim_component(coordinate_index))
                .unwrap()
        });
        let keys = nodes.iter().map(|n| &n.centroid).collect::<Vec<_>>();

        let mut min_variance = f32::INFINITY;
        let mut split_index = m;
        for i in m..=(nodes.len() - m) {
            let variance1 = K::variance_along_direction(&keys[..i], coordinate_index);
            let variance2 = K::variance_along_direction(&keys[i..], coordinate_index);
            let variance = variance1 + variance2;
            if variance < min_variance {
                min_variance = variance;
                split_index = i;
            }
        }
        split_index
    }

    pub fn find_sibling_to_borrow_from<
        P,
        K: Distance + MaxVarianceDirection + PartialEq,
        const M: usize,
    >(
        nodes: &[SsNode<P, K, M>],
        node_to_fix: usize,
        m: usize,
    ) -> Option<usize> {
        let siblings_to_borrow_from =
            nodes
                .iter()
                .enumerate()
                .filter(|(i, sibling)| match &sibling.links {
                    SsNodeLinks::Inner(nodes) => *i != node_to_fix && nodes.len() > m,
                    SsNodeLinks::Leaf(points) => *i != node_to_fix && points.len() > m,
                });

        let mut closest_sibling = None;
        let mut closest_sibling_dist = f32::INFINITY;

        for (i, sibling) in siblings_to_borrow_from {
            let distance = nodes[node_to_fix].centroid.distance(&sibling.centroid);
            if distance < closest_sibling_dist {
                closest_sibling = Some(i);
                closest_sibling_dist = distance;
            }
        }
        closest_sibling
    }

    pub fn borrow_from_sibling<
        P,
        K: Distance + Zero + MaxVarianceDirection + PartialEq,
        const M: usize,
    >(
        nodes: &mut [SsNode<P, K, M>],
        node_to_fix: usize,

        sibling_to_borrow_from: usize,
    ) {
        // found sibling to borrow from
        let to_fix_centroid = &nodes[node_to_fix].centroid;

        match &mut nodes[sibling_to_borrow_from].links {
            SsNodeLinks::Inner(nodes2) => {
                let mut closest_node = None;
                let mut closest_node_dist = f32::INFINITY;
                for (i, node) in nodes2.iter().enumerate() {
                    let distance = node.centroid.distance(to_fix_centroid);
                    if distance < closest_node_dist {
                        closest_node = Some(i);
                        closest_node_dist = distance;
                    }
                }
                let node = nodes2.remove(closest_node.unwrap());
                nodes[sibling_to_borrow_from].update_bounding_envelope();

                match &mut nodes[node_to_fix].links {
                    SsNodeLinks::Inner(fix_nodes) => fix_nodes.push(node),
                    SsNodeLinks::Leaf(_) => panic!("unbalanced tree"),
                }
                nodes[node_to_fix].update_bounding_envelope();
            }
            SsNodeLinks::Leaf(points) => {
                let mut closest_point = None;
                let mut closest_point_dist = f32::INFINITY;
                for (i, point) in points.iter().enumerate() {
                    let distance = point.center.distance(to_fix_centroid);
                    if distance < closest_point_dist {
                        closest_point = Some(i);
                        closest_point_dist = distance;
                    }
                }
                // println!(
                //     "closest point: {:?} {} {}",
                //     closest_point, sibling_to_borrow_from, node_to_fix
                // );
                let point = points.remove(closest_point.unwrap());
                nodes[sibling_to_borrow_from].update_bounding_envelope();
                match &mut nodes[node_to_fix].links {
                    SsNodeLinks::Inner(_) => panic!("unbalanced tree"),
                    SsNodeLinks::Leaf(fix_points) => fix_points.push(point),
                }
                nodes[node_to_fix].update_bounding_envelope();
            }
        }
    }

    pub fn find_sibling_to_merge_to<
        P,
        K: Distance + MaxVarianceDirection + PartialEq,
        const M: usize,
    >(
        nodes: &[SsNode<P, K, M>],
        node_to_fix: usize,
        m: usize,
    ) -> Option<usize> {
        let siblings_to_merge_to =
            nodes
                .iter()
                .enumerate()
                .filter(|(i, sibling)| match &sibling.links {
                    SsNodeLinks::Inner(nodes) => *i != node_to_fix && nodes.len() == m,
                    SsNodeLinks::Leaf(points) => *i != node_to_fix && points.len() == m,
                });

        let mut closest_sibling = None;
        let mut closest_sibling_dist = f32::INFINITY;

        for (i, sibling) in siblings_to_merge_to {
            let distance = nodes[node_to_fix].centroid.distance(&sibling.centroid);
            if distance < closest_sibling_dist {
                closest_sibling = Some(i);
                closest_sibling_dist = distance;
            }
        }
        closest_sibling
    }

    pub fn merge_siblings<
        P,
        K: Zero + Distance + MaxVarianceDirection + PartialEq,
        const M: usize,
    >(
        nodes: &mut ArrayVec<SsNode<P, K, M>, M>,
        mut node_index_1: usize,
        mut node_index_2: usize,
    ) {
        if node_index_1 > node_index_2 {
            // remove node with larger index first
            std::mem::swap(&mut node_index_1, &mut node_index_2);
        }
        let node_2 = nodes.remove(node_index_2);
        let node_1 = nodes.remove(node_index_1);
        let node = merge(node_1, node_2);
        nodes.push(node);
    }

    fn merge<P, K: Zero + Distance + MaxVarianceDirection + PartialEq, const M: usize>(
        node_1: SsNode<P, K, M>,
        node_2: SsNode<P, K, M>,
    ) -> SsNode<P, K, M> {
        match (node_1.links, node_2.links) {
            (SsNodeLinks::Leaf(mut points1), SsNodeLinks::Leaf(mut points2)) => {
                points1.extend(points2.drain(..));
                SsNode::<P, K, M>::from_elements(*points1)
            }
            (SsNodeLinks::Inner(mut nodes1), SsNodeLinks::Inner(mut nodes2)) => {
                nodes1.extend(nodes2.drain(..));
                SsNode::<P, K, M>::from_nodes(*nodes1)
            }
            _ => panic!("inconsistent siblings"),
        }
    }
}
#[test]
fn test_search() {
    const UPPER_M: usize = 8;
    const LOWER_M: usize = 4;

    let mut tree = SsTree::<(), [f32; 2], UPPER_M>::new(LOWER_M);

    tree.insert(Element::new([0.0, 0.0], 1.0, ()));
    tree.insert(Element::new([5.0, 5.0], 1.0, ()));

    let mut out = Vec::new();
    tree.points_within_radius(&[0.5, 0.5], 1.0, &mut out);
    assert_eq!(out, vec!(&Element::new([0.0, 0.0], 1.0, ())));

    let mut out = Vec::new();
    tree.points_within_radius(&[4.5, 5.5], 1.0, &mut out);
    assert_eq!(out, vec!(&Element::new([5.0, 5.0], 1.0, ())));
    let mut out = Vec::new();

    // do search between the elements with radius big enough to just reach them
    tree.points_within_radius(
        &[2.5, 2.5],
        (2.5 * std::f32::consts::SQRT_2 + 0.0001) - 1.0,
        &mut out,
    );
    assert_eq!(out.len(), 2);
    assert!(out.contains(&&Element::new([5.0, 5.0], 1.0, ())));
    assert!(out.contains(&&Element::new([0.0, 0.0], 1.0, ())));

    let mut out = Vec::new();

    // the same as befor but with radius just barely too small
    tree.points_within_radius(
        &[2.5, 2.5],
        (2.5 * std::f32::consts::SQRT_2 - 0.0001) - 1.0,
        &mut out,
    );
    assert!(out.is_empty());
}
