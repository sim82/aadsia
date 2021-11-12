use std::ops::Index;

use arrayvec::ArrayVec;

#[derive(Debug)]
pub enum SsNodeLinks<const K: usize, const M: usize> {
    Inner(ArrayVec<Box<SsNode<K, M>>, M>),
    Leaf(ArrayVec<[f32; K], M>),
}

#[derive(Debug)]
pub struct SsNode<const K: usize, const M: usize> {
    pub centroid: [f32; K],
    pub radius: f32,
    pub links: SsNodeLinks<K, M>,
}

fn distance<const K: usize>(p1: &[f32; K], p2: &[f32; K]) -> f32 {
    p1.iter()
        .zip(p2.iter())
        .map(|(c1, c2)| (c1 - c2) * (c1 - c2))
        .sum::<f32>()
        .sqrt()
}

#[test]
fn test_distance() {
    assert_eq!(distance(&[0.0, 0.0], &[1.0, 1.0]), 2.0f32.sqrt());
    assert_eq!(distance(&[-10.0, 1.0], &[10.0, 1.0]), 20.0f32);
    assert_eq!(distance(&[1000.0, -1000.0], &[1000.0, 2000.0]), 3000.0);
}

impl<const K: usize, const M: usize> SsNode<K, M> {
    pub fn intersects_point(&self, target: &[f32; K]) -> bool {
        distance(&self.centroid, target) <= self.radius
    }

    pub fn search(&self, target: &[f32; K]) -> Option<&Self> {
        match &self.links {
            SsNodeLinks::Inner(children) => {
                children
                    .iter()
                    .find(|node| node.intersects_point(target))
                    .map(|node| node.as_ref())
                // for node in children {
                //     if node.intersects_points(target) {
                //         match node.search(target) {
                //             Some(node) => return Some(node),
                //             None => continue,
                //         }
                //     }
                // }
                // return None;
            }
            SsNodeLinks::Leaf(points) => {
                if points.iter().any(|x| *x == *target) {
                    Some(self)
                } else {
                    None
                }
            }
            _ => panic!("split link types not allowed"),
        }
    }

    pub fn search_parent_leaf(&self, target: &[f32; K]) -> &Self {
        match &self.links {
            SsNodeLinks::Inner(children) => {
                let child = find_closest_child(children, target);
                child.search_parent_leaf(target)
            }
            SsNodeLinks::Leaf(_) => self,
            _ => panic!("split link types not allowed"),
        }
    }

    pub fn update_bounding_envelope(&mut self) {
        let (centroid, radius) = match &self.links {
            SsNodeLinks::Inner(nodes) => inner::centroid_and_radius(nodes),
            SsNodeLinks::Leaf(points) => leaf::centroid_and_radius(points),
        };
        self.centroid = centroid;
        self.radius = radius;
    }
    pub fn insert(&mut self, point: &[f32; K], m: usize) -> Option<(Box<Self>, Box<Self>)> {
        match &mut self.links {
            SsNodeLinks::Leaf(points) => {
                if points.iter().any(|p| *p == *point) {
                    return None;
                }

                if points.len() < M {
                    points.push(*point);
                    self.update_bounding_envelope();
                    return None;
                } else {
                    let mut nodes_to_split: Vec<[f32; K]> =
                        points.drain(..).chain(std::iter::once(*point)).collect();

                    let split_index = leaf::find_split_index(&mut nodes_to_split, m);
                    let points1: ArrayVec<_, M> =
                        nodes_to_split[..split_index].iter().cloned().collect();
                    let (centroid1, radius1) = leaf::centroid_and_radius(&points1);

                    let points2: ArrayVec<_, M> =
                        nodes_to_split[split_index..].iter().cloned().collect();
                    let (centroid2, radius2) = leaf::centroid_and_radius(&points2);

                    let new_node1 = Box::new(Self {
                        centroid: centroid1,
                        radius: radius1,
                        links: SsNodeLinks::Leaf(points1),
                    });
                    let new_node2 = Box::new(Self {
                        centroid: centroid2,
                        radius: radius2,
                        links: SsNodeLinks::Leaf(points2),
                    });
                    return Some((new_node1, new_node2));
                }
            }

            SsNodeLinks::Inner(children) => {
                let closest_child_index = find_closest_child_index(children, point);
                if let Some((new_child_1, new_child_2)) =
                    children[closest_child_index].insert(point, m)
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

                        let new_node1 = Box::new(Self {
                            centroid: centroid1,
                            radius: radius1,
                            links: SsNodeLinks::Inner(points1),
                        });
                        let new_node2 = Box::new(Self {
                            centroid: centroid2,
                            radius: radius2,
                            links: SsNodeLinks::Inner(points2),
                        });
                        return Some((new_node1, new_node2));
                    }
                } else {
                    self.update_bounding_envelope();
                }
            }
            _ => panic!("split link types not allowed"),
        }
        None
        // self.split()
        //   if this.leaf then
        //     if point in this.points then
        //       return null
        //     this.points.add(point)
        //     this.updateBoundingEnvelope()
        //     if this.points.size <= M then
        //       return null
        //   else
        //     closestChild ← this.findClosestChild()
        //     (newChild1, newChild2) ← insert(closestChild, point)
        //     if newChild1 == null then
        //       node.updateBoundingEnvelope()
        //       return null
        //     else
        //       this.children.delete(closestChild)
        //       this.children.add(newChild1)
        //       this.children.add(newChild2)
        //       node.updateBoundingEnvelope()
        //       if this.children.size <= M then
        //         return null
        //   return this.split()
    }

    pub fn delete(&mut self, target: &[f32; K], m: usize) -> (bool, bool) {
        match &mut self.links {
            SsNodeLinks::Leaf(points) => {
                if let Some((i, _)) = points.iter().enumerate().find(|(_, p)| *p == target) {
                    points.remove(i);
                    (true, points.len() < m)
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
                        let siblings_to_borrow_from =
                            nodes
                                .iter()
                                .enumerate()
                                .filter(|(i, sibling)| match &sibling.links {
                                    SsNodeLinks::Inner(nodes) => {
                                        *i != node_to_fix && nodes.len() > m
                                    }
                                    SsNodeLinks::Leaf(points) => {
                                        *i != node_to_fix && points.len() > m
                                    }
                                });

                        let mut closest_sibling = None;
                        let mut closest_sibling_dist = f32::INFINITY;

                        for (i, sibling) in siblings_to_borrow_from {
                            let distance =
                                distance(&nodes[node_to_fix].centroid, &sibling.centroid);
                            if distance < closest_sibling_dist {
                                closest_sibling = Some(i);
                                closest_sibling_dist = distance;
                            }
                        }

                        if let Some(closest_sibling) = closest_sibling {
                            let to_fix_centroid = nodes[node_to_fix].centroid;

                            match &mut nodes[closest_sibling].links {
                                SsNodeLinks::Inner(nodes2) => {
                                    let mut closest_node = None;
                                    let mut closest_node_dist = f32::INFINITY;
                                    for (i, node) in nodes2.iter().enumerate() {
                                        let distance = distance(&node.centroid, &to_fix_centroid);
                                        if distance < closest_node_dist {
                                            closest_node = Some(i);
                                            closest_node_dist = distance;
                                        }
                                    }
                                    let node = nodes2.remove(closest_node.unwrap());
                                    nodes[closest_sibling].update_bounding_envelope();

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
                                        let distance = distance(point, &to_fix_centroid);
                                        if distance < closest_point_dist {
                                            closest_point = Some(i);
                                            closest_point_dist = distance;
                                        }
                                    }
                                    println!(
                                        "closest point: {:?} {} {}",
                                        closest_point, closest_sibling, node_to_fix
                                    );
                                    let point = points.remove(closest_point.unwrap());
                                    nodes[closest_sibling].update_bounding_envelope();
                                    match &mut nodes[node_to_fix].links {
                                        SsNodeLinks::Inner(_) => panic!("unbalanced tree"),
                                        SsNodeLinks::Leaf(fix_points) => fix_points.push(point),
                                    }
                                    nodes[node_to_fix].update_bounding_envelope();
                                }
                            }
                        }

                        (true, nodes.len() < m)
                    }
                }
            }
        }
    }

    //     function delete(node, target)
    //   if node.leaf then
    //     if node.points.contains(target) then
    //       node.points.delete(target)
    //       return (true, node.points.size() < m)
    //     else
    //      return (false, false)
    //   else
    //     nodeToFix ← null
    //     deleted ← false
    //     for childNode in node.children do
    //       if childNode.intersectsPoint(target) then
    //         (deleted, violatesevariants) ← delete(childNode, target)
    //         if violatesevariants == true then
    //           nodeToFix ← childNode
    //         if deleted then
    //           break
    //   if nodeToFix == null then
    //     if deleted then
    //       node.updateBoundingEnvelope()
    //     return (deleted, false)
    //   else
    //     siblingsToBorrowFrom(nodeToFix)
    //     if not siblings.isEmpty() then
    //       nodeToFix.borrowFromSibling(siblings)
    //     else
    //       node.mergeChildren(
    //         nodeToFix, node.findSiblingToMergeTo(nodeToFix))
    //     node.updateBoundingEnvelope()
    //     return (true, node.children.size() < m)

    pub fn count_nodes(&self) -> (usize, usize) {
        match &self.links {
            SsNodeLinks::Inner(nodes) => nodes.iter().fold((0, 1), |(a_points, a_nodes), n| {
                let (points, nodes) = n.count_nodes();
                (a_points + points, a_nodes + nodes)
            }),
            SsNodeLinks::Leaf(points) => (points.len(), 1),
        }
    }
}
fn find_closest_child<'a, const K: usize, const M: usize>(
    children: &'a [Box<SsNode<K, M>>],
    target: &[f32; K],
) -> &'a SsNode<K, M> {
    // children
    //     .iter()
    //     .min_by_key(|a| distance(&a.centroid, target))
    //     .as_ref()
    //     .unwrap()

    let mut min_dist = f32::MAX;
    let mut cur_min = None;
    for child in children {
        let d = distance(&child.centroid, target);
        if d < min_dist {
            min_dist = d;
            cur_min = Some(child.as_ref());
        }
    }
    cur_min.unwrap()
}
fn find_closest_child_index<const K: usize, const M: usize>(
    children: &[Box<SsNode<K, M>>],
    target: &[f32; K],
) -> usize {
    let mut min_dist = f32::MAX;
    let mut cur_min = None;
    for (i, child) in children.iter().enumerate() {
        let d = distance(&child.centroid, target);
        if d < min_dist {
            min_dist = d;
            cur_min = Some(i);
        }
    }
    cur_min.unwrap()
}

#[derive(Debug)]
pub struct SsTree<const K: usize, const M: usize> {
    pub root: SsNode<K, M>,
    height: usize,
    m: usize,
}

impl<const K: usize, const M: usize> SsTree<K, M> {
    pub fn new(m: usize) -> Self {
        Self {
            root: SsNode {
                centroid: [0f32; K],
                radius: 0f32,
                links: SsNodeLinks::Leaf(ArrayVec::new()),
            },
            height: 1,
            m,
        }
    }

    pub fn insert(&mut self, point: &[f32; K]) {
        if let Some((new_child_1, new_child_2)) = self.root.insert(point, self.m) {
            let mut nodes = ArrayVec::<_, M>::new();
            nodes.push(new_child_1);
            nodes.push(new_child_2);
            let (centroid, radius) = inner::centroid_and_radius(&nodes);
            self.root = SsNode {
                centroid,
                radius,
                links: SsNodeLinks::Inner(nodes),
            };
            self.height += 1;
        }
    }
    pub fn delete(&mut self, point: &[f32; K]) {
        self.root.delete(point, self.m);
    }

    pub fn get_height(&self) -> usize {
        self.height
    }
    pub fn get_fill_factor(&self) -> f32 {
        let (num_points, num_nodes) = self.root.count_nodes();
        num_points as f32 / num_nodes as f32
    }
}

impl<const K: usize, const M: usize> Default for SsTree<K, M> {
    fn default() -> Self {
        Self::new(M / 2)
    }
}

mod leaf {
    pub fn mean_along_direction<const K: usize>(
        points: &[[f32; K]],
        direction_index: usize,
    ) -> f32 {
        assert!(!points.is_empty());
        let count = points.len() as f32;
        let sum = points
            .iter()
            .map(|point| point[direction_index])
            .sum::<f32>();
        sum / count
    }

    pub fn variance_along_direction<const K: usize>(
        points: &[[f32; K]],
        direction_index: usize,
    ) -> f32 {
        assert!(!points.is_empty());
        let mean = mean_along_direction(points, direction_index);
        let count = points.len() as f32;
        points
            .iter()
            .map(|point| {
                let diff = mean - point[direction_index];

                diff * diff
            })
            .sum::<f32>()
            / count
    }

    pub fn direction_of_max_variance<const K: usize>(points: &[[f32; K]]) -> usize {
        let mut max_variance = 0.0;
        let mut direction_index = 0;
        for i in 0..K {
            let variance = variance_along_direction(points, i);
            if variance > max_variance {
                max_variance = variance;
                direction_index = i;
            }
        }
        direction_index
    }

    pub fn find_split_index<const K: usize>(points: &mut [[f32; K]], m: usize) -> usize {
        let coordinate_index = direction_of_max_variance(points);
        points.sort_by(|p1, p2| {
            p1[coordinate_index]
                .partial_cmp(&p2[coordinate_index])
                .unwrap()
        });
        let mut min_variance = f32::INFINITY;
        let mut split_index = m;
        for i in m..=(points.len() - m) {
            let variance1 = variance_along_direction(&points[..i], coordinate_index);
            let variance2 = variance_along_direction(&points[i..], coordinate_index);
            let variance = variance1 + variance2;
            if variance < min_variance {
                min_variance = variance;
                split_index = i;
            }
        }
        split_index
    }

    pub fn centroid_and_radius<const K: usize>(points: &[[f32; K]]) -> ([f32; K], f32) {
        let mut centroid = [0f32; K];
        for i in 0..K {
            centroid[i] = mean_along_direction(points, i);
        }

        let radius = points
            .iter()
            .map(|point| super::distance(&centroid, point))
            .max_by(|d1, d2| d1.partial_cmp(d2).unwrap())
            .unwrap();
        (centroid, radius)
    }
}
mod inner {
    use super::{distance, SsNode};

    pub fn find_split_index<const K: usize, const M: usize>(
        nodes: &mut [Box<SsNode<K, M>>],
        m: usize,
    ) -> usize {
        let coordinate_index = direction_of_max_variance(nodes);
        nodes.sort_by(|p1, p2| {
            p1.centroid[coordinate_index]
                .partial_cmp(&p2.centroid[coordinate_index])
                .unwrap()
        });
        let mut min_variance = f32::INFINITY;
        let mut split_index = m;
        for i in m..=(nodes.len() - m) {
            let variance1 = variance_along_direction(&nodes[..i], coordinate_index);
            let variance2 = variance_along_direction(&nodes[i..], coordinate_index);
            let variance = variance1 + variance2;
            if variance < min_variance {
                min_variance = variance;
                split_index = i;
            }
        }
        split_index
    }

    pub fn centroid_and_radius<const K: usize, const M: usize>(
        nodes: &[Box<SsNode<K, M>>],
    ) -> ([f32; K], f32) {
        let mut centroid = [0f32; K];
        for i in 0..K {
            centroid[i] = mean_along_direction(nodes, i);
        }

        // let centroid = mean_along_all_directions(nodes);

        let radius = nodes
            .iter()
            .map(|node| distance(&centroid, &node.centroid) + node.radius)
            .max_by(|d1, d2| d1.partial_cmp(d2).unwrap())
            .unwrap();
        (centroid, radius)
    }

    pub fn mean_along_direction<const K: usize, const M: usize>(
        points: &[Box<SsNode<K, M>>],
        direction_index: usize,
    ) -> f32 {
        assert!(!points.is_empty());
        let count = points.len() as f32;
        let sum = points
            .iter()
            .map(|point| point.centroid[direction_index])
            .sum::<f32>();
        sum / count
    }

    // pub fn mean_along_all_directions<const K: usize, const M: usize>(
    //     points: &[Box<SsNode<K, M>>],
    // ) -> [f32; K] {
    //     assert!(!points.is_empty());
    //     let count = points.len() as f32;
    //     let mut sum = points.iter().fold([0f32; K], |a, node| {
    //         let mut a = a.clone();
    //         a.iter_mut()
    //             .zip(node.centroid.iter())
    //             .for_each(|(a, c)| *a += *c);
    //         a
    //     });
    //     sum.iter_mut().for_each(|sum| *sum /= count);
    //     sum
    // }

    pub fn variance_along_direction<const K: usize, const M: usize>(
        points: &[Box<SsNode<K, M>>],
        direction_index: usize,
    ) -> f32 {
        assert!(!points.is_empty());
        let mean = mean_along_direction(points, direction_index);
        let count = points.len() as f32;
        points
            .iter()
            .map(|point| {
                let diff = mean - point.centroid[direction_index];

                diff * diff
            })
            .sum::<f32>()
            / count
    }

    pub fn direction_of_max_variance<const K: usize, const M: usize>(
        points: &[Box<SsNode<K, M>>],
    ) -> usize {
        let mut max_variance = 0.0;
        let mut direction_index = 0;
        for i in 0..K {
            let variance = variance_along_direction(points, i);
            if variance > max_variance {
                max_variance = variance;
                direction_index = i;
            }
        }
        direction_index
    }
}
#[test]
fn test_search() {
    // let root =
}
