use arrayvec::ArrayVec;

enum SsNodeLinks<const K: usize, const M: usize> {
    Inner(ArrayVec<Box<SsNode<K, M>>, M>),
    Leaf(ArrayVec<[f32; K], M>),
    SplitInner(Vec<Box<SsNode<K, M>>>),
    SplitLeaf(Vec<[f32; K]>),
}

struct SsNode<const K: usize, const M: usize> {
    centroid: [f32; K],
    radius: f32,
    links: SsNodeLinks<K, M>,
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
    pub fn intersects_points(&self, target: &[f32; K]) -> bool {
        distance(&self.centroid, target) <= self.radius
    }

    pub fn search(&self, target: &[f32; K]) -> Option<&Self> {
        match &self.links {
            SsNodeLinks::Inner(children) => {
                children
                    .iter()
                    .find(|node| node.intersects_points(target))
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

    pub fn update_bounding_envelope(&mut self) -> Option<&Self> {
        None
    }
    pub fn insert(&mut self, point: &[f32; K]) -> Option<&Self> {
        match &mut self.links {
            SsNodeLinks::Leaf(points) => {
                if points.iter().any(|p| *p == *point) {
                    return None;
                }

                if points.len() < M - 2 {
                    points.push(*point);
                    self.update_bounding_envelope();
                    return None;
                } else {
                    self.links = SsNodeLinks::SplitLeaf(
                        points.drain(..).chain(std::iter::once(*point)).collect(),
                    )
                }
            }

            SsNodeLinks::Inner(_) => todo!(),
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

#[test]
fn test_search() {
    // let root =
}
