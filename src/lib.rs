#[derive(Debug, PartialEq)]
pub struct KdNode<const N: usize> {
    point: [f32; N],
    left: Option<Box<KdNode<N>>>,
    right: Option<Box<KdNode<N>>>,
    level: usize,
}

pub fn get_point_key<const N: usize>(point: &[f32; N], level: usize) -> f32 {
    point[level % N]
}

impl<const N: usize> KdNode<N> {
    pub fn new(
        point: [f32; N],
        left: Option<KdNode<N>>,
        right: Option<KdNode<N>>,
        level: usize,
    ) -> Self {
        KdNode {
            point,
            left: left.map(|left| Box::new(left)),
            right: right.map(|right| Box::new(right)),
            level,
        }
    }

    pub fn get_key(&self) -> f32 {
        get_point_key(&self.point, self.level)
    }
    pub fn compare(&self, point: &[f32; N]) -> f32 {
        (get_point_key(point, self.level) - self.get_key()).signum()
    }
    pub fn split_distance(&self, point: &[f32; N]) -> f32 {
        (get_point_key(point, self.level) - self.get_key()).abs()
    }
    pub fn search(&self, target: &[f32; N]) -> Option<&KdNode<N>> {
        self.dump();
        if self.point == *target {
            Some(self)
        } else if self.compare(target) < 0.0 {
            self.left.as_ref()?.search(target)
        } else {
            self.right.as_ref()?.search(target)
        }
    }
    pub fn dump(&self) {
        println!(
            "{}level: {}, point:{:?}",
            std::iter::repeat(' ').take(self.level).collect::<String>(),
            self.level,
            self.point
        )
    }
    pub fn dump_rec(&self, p: &str) {
        println!(
            "{}{}level: {}, point:{:?}",
            std::iter::repeat(' ').take(self.level).collect::<String>(),
            p,
            self.level,
            self.point
        );

        if let Some(ref left) = self.left {
            left.dump_rec("l ");
        }
        if let Some(ref right) = self.right {
            right.dump_rec("r ");
        }
    }

    pub fn insert(&mut self, new_point: &[f32; N], level: usize) -> &KdNode<N> {
        if self.point == *new_point {
            self
        } else if self.compare(new_point) < 0.0 {
            match self.left.as_mut() {
                Some(left) => {
                    left.insert(new_point, level + 1);
                }
                None => {
                    self.left = Some(Box::new(KdNode::<N>::new(
                        *new_point,
                        None,
                        None,
                        level + 1,
                    )))
                }
            };
            return self;
            // match node.left {
            //     Some(left) => left.insert(new_point, level + 1),
            //     None => node.left.
            // }
        } else {
            match self.right.as_mut() {
                Some(right) => {
                    right.insert(new_point, level + 1);
                }
                None => {
                    self.right = Some(Box::new(KdNode::<N>::new(
                        *new_point,
                        None,
                        None,
                        level + 1,
                    )))
                }
            };
            return self;
        }
        // if node == null then
        //   return new KdNode(newPoint, null, null, level)
        // elsif node.point == newPoint then
        //   return node
        // elsif compare(newPoint, node) < 0 then
        //   node.left ← insert(node.left, newPoint, node.level + 1)
        //   return node
        // else
        //   node.right ← insert(node.right, newPoint, node.level + 1)
        //   return node
    }
}

pub fn partition<const N: usize>(
    points: &mut [[f32; N]],
    level: usize,
) -> ([f32; N], &mut [[f32; N]], &mut [[f32; N]]) {
    assert!(!points.is_empty());
    if points.len() == 1 {
        return (points[0], Default::default(), Default::default());
    }
    points.sort_by(|point1, point2| {
        get_point_key(point1, level)
            .partial_cmp(&get_point_key(point2, level))
            .unwrap()
    });

    if points.len() == 2 {
        return (points[1], &mut points[0..0], Default::default());
    }
    let m = points.len() / 2;
    let (left, right) = points.split_at_mut(m);
    let (m, right) = right.split_at_mut(1);
    (m[0], left, right)
}

pub fn construct_balanced<const N: usize>(
    points: &mut [[f32; N]],
    level: usize,
) -> Option<KdNode<N>> {
    if points.is_empty() {
        None
    } else if points.len() == 1 {
        Some(KdNode::new(points[0], None, None, level))
    } else {
        let (median, left, right) = partition(points, level);
        let left_tree = construct_balanced(left, level + 1);
        let right_tree = construct_balanced(right, level + 1);
        Some(KdNode::new(median, left_tree, right_tree, level))
    }
}

struct KdTree<const N: usize> {
    root: KdNode<N>,
    k: usize,
}

impl<const N: usize> KdTree<N> {
    // pub fn new(points: &[[f32; N]]) -> Self {}
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        type KdNode2 = KdNode<2>;

        let mut root = KdNode2::new(
            [0.0, 5.0],
            Some(KdNode2::new(
                [-1.0, 6.0],
                Some(KdNode2::new(
                    [-1.0, 1.0],
                    None,
                    Some(KdNode2::new([-0.5, 0.0], None, None, 3)),
                    2,
                )),
                None,
                1,
            )),
            Some(KdNode2::new(
                [1.0, -1.0],
                Some(KdNode2::new([2.0, -5.0], None, None, 2)),
                None,
                1,
            )),
            0,
        );

        let target1 = KdNode2::new(
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
                Some(KdNode2::new([2.0, -5.0], None, None, 2)),
                None,
                1,
            )),
            0,
        );

        let target2 = KdNode2::new(
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

        root.insert(&[-1.5, -2.0], 0);
        assert_eq!(root, target1);
        root.insert(&[2.5, -3.0], 0);
        assert_eq!(root, target2);
    }

    #[test]
    fn test_partition() {
        let mut points = [
            [1.0, -1.0],
            [-1.0, 6.0],
            [-0.5, 0.0],
            [2.0, 5.0],
            [0.0, 5.0],
            [2.5, 3.0],
            [-1.0, 1.0],
            [-1.5, -2.0],
        ];

        let (median, left, right) = partition(&mut points, 0);
        println!("{:?} {:?} {:?}", median, left, right);
    }

    #[test]
    fn test_construct_balanced() {
        let mut points = [
            [1.0, -1.0],
            [-1.0, 6.0],
            [-0.5, 0.0],
            [2.0, 5.0],
            [0.0, 5.0],
            [2.5, 3.0],
            [-1.0, 1.0],
            [-1.5, -2.0],
        ];

        let tree = construct_balanced(&mut points, 0);
        tree.unwrap().dump_rec("");
        // println!("{:?}", tree);
    }
}
