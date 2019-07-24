// Copyright (c) 2015, Mikhail Vorotilov
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use super::super::FloatType;
use core::iter::FusedIterator;

#[derive(Default)]
pub struct Roots<F> where F: FloatType {
    roots: [F; 4],
    num_roots: usize,
    cursor: usize,

}

impl<F> Roots<F> where F: FloatType {
    pub fn add_new_root(&mut self, root: F) {
        if self.num_roots < self.roots.len() {
            // let i = {
            //     let mut i = 0;
            //     while i < self.num_roots {
            //         if root < self.roots[i] {
            //             break;
            //         }

            //         i += 1;
            //     }
            //     i
            // };

            // unsafe { 
            //     self.roots[i..]
            //         .as_mut_ptr()
            //         .copy_to(self.roots[i+1..].as_mut_ptr(), self.num_roots - i);
            // }
            
            // self.roots[i] = root;
            // self.num_roots += 1;
        }
    }   

    pub fn zero() -> Self {
        Self::default()
    }

    pub fn one(root: F) -> Self {
        Self {
            roots: [root, F::default(), F::default(), F::default()],
            num_roots: 1,
            ..Default::default()
        }
    }

    pub fn two(root1: F, root2: F) -> Self {
        Self {
            roots: [root1, root2, F::default(), F::default()],
            num_roots: 2,
            ..Default::default()
        }
    }

    pub fn three(root1: F, root2: F, root3: F) -> Self {
        Self {
            roots: [root1, root2, root3, F::default()],
            num_roots: 3,
            ..Default::default()
        }
    }

    pub fn four(root1: F, root2: F, root3: F, root4: F) -> Self {
        Self {
            roots: [root1, root2, root3, root4],
            num_roots: 4,
            ..Default::default()
        }
    }
}

impl<F> Iterator for Roots<F> where F: FloatType {
    type Item = F;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor == self.num_roots {
            None
        } else {
            let root = self.roots[self.cursor];
            self.cursor += 1;
            Some(root)
        }
    }
}

impl<F> FusedIterator for Roots<F> where F: FloatType {}