#[derive(Debug, PartialEq)]
struct KdNode<const N: usize> {
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
        let mut root = KdNode::<2> {
            point: [0.0, 5.0],
            left: Some(Box::new(KdNode::<2> {
                point: [-1.0, 6.0],
                left: Some(Box::new(KdNode::<2> {
                    point: [-1.0, 1.0],
                    left: None,
                    right: Some(Box::new(KdNode::<2> {
                        point: [-0.5, 0.0],
                        left: None,
                        right: None,
                        level: 3,
                    })),
                    level: 2,
                })),
                right: None,
                level: 1,
            })),
            right: Some(Box::new(KdNode::<2> {
                point: [1.0, -1.0],
                left: Some(Box::new(KdNode::<2> {
                    point: [2.0, -5.0],
                    left: None,
                    right: None,
                    level: 2,
                })),
                right: None,
                level: 1,
            })),
            level: 0,
        };

        let target1 = KdNode::<2> {
            point: [0.0, 5.0],
            left: Some(Box::new(KdNode::<2> {
                point: [-1.0, 6.0],
                left: Some(Box::new(KdNode::<2> {
                    point: [-1.0, 1.0],
                    left: Some(Box::new(KdNode::<2> {
                        point: [-1.5, -2.0],
                        left: None,
                        right: None,
                        level: 3,
                    })),
                    right: Some(Box::new(KdNode::<2> {
                        point: [-0.5, 0.0],
                        left: None,
                        right: None,
                        level: 3,
                    })),
                    level: 2,
                })),
                right: None,
                level: 1,
            })),
            right: Some(Box::new(KdNode::<2> {
                point: [1.0, -1.0],
                left: Some(Box::new(KdNode::<2> {
                    point: [2.0, -5.0],
                    left: None,
                    right: None,
                    level: 2,
                })),
                right: None,
                level: 1,
            })),
            level: 0,
        };

        let target2 = KdNode::<2> {
            point: [0.0, 5.0],
            left: Some(Box::new(KdNode::<2> {
                point: [-1.0, 6.0],
                left: Some(Box::new(KdNode::<2> {
                    point: [-1.0, 1.0],
                    left: Some(Box::new(KdNode::<2> {
                        point: [-1.5, -2.0],
                        left: None,
                        right: None,
                        level: 3,
                    })),
                    right: Some(Box::new(KdNode::<2> {
                        point: [-0.5, 0.0],
                        left: None,
                        right: None,
                        level: 3,
                    })),
                    level: 2,
                })),
                right: None,
                level: 1,
            })),
            right: Some(Box::new(KdNode::<2> {
                point: [1.0, -1.0],
                left: Some(Box::new(KdNode::<2> {
                    point: [2.0, -5.0],
                    left: None,
                    right: Some(Box::new(KdNode::<2> {
                        point: [2.5, -3.0],
                        left: None,
                        right: None,
                        level: 3,
                    })),
                    level: 2,
                })),
                right: None,
                level: 1,
            })),
            level: 0,
        };

        root.insert(&[-1.5, -2.0], 0);
        assert_eq!(root, target1);
        root.insert(&[2.5, -3.0], 0);
        assert_eq!(root, target2);
    }
}
