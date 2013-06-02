-- | A specification for a Support Vector Machine algorithm. VERY SLOW!
{-# OPTIONS_GHC -fwarn-missing-signatures #-}

-- Copyright Ezra Cooper, circa 2008.

{-
  TODO:
    * Choose the candidates more cleverly.
    * Use Data.Vector for matrix representations.  
    * Examine the Gaussian algorithm. Faster versions?
    * Performance benchmarking.
 -}

module SVM (margin, classify, example, example2, example3) where

import Control.Monad (MonadPlus, mplus, mzero, msum, replicateM)
import Data.List (tails, partition)
import System.Random
import Test.HUnit

import Debug.Trace

-- Local imports:
import Matrices

-- | The type of linear classifiers.
data LC a  =  LC(Vector a -- ^ weight vector.
                , a       -- ^ cutoff point along that vector.
                )
    deriving Show

-- | Given a (binary) linear classifier @l@ and an example vector @v@,
-- which class @l@ give to @v@?
classify :: (Ord a, Floating a) => LC a -> Vector a -> Bool
classify (LC(w, b)) x = w % x - b <= 0

-- | Find the maximal-margin hyperplane between two datasets (tagged as True and False in the given list).
margin :: [(Bool, Vector Double)] -> LC Double
margin [] = error "Need some bloody points!"
margin pts = 
  let dim = 1+length (snd $ head pts) in
  if length pts < dim then error "Too few points for unambiguous separation"
  else
  -- | The adjusted points: tag all points with {1, -1} and invert one of the
  -- sets in space.
  let adjPts = map (\(c, x) -> if c then 1 : neg x else (-1) : x) pts in
  -- | Hyperplanes that might define the boundary, because they cut through
  -- @dim@ points from the adjusted set.
  let candidates = joinMaybes [gaussianElimUniq (plane <|> repeat 1)
                               | plane <- combos dim adjPts] in
  -- | Test that all adjPts are "beyond" the given vector's hyperplane.
  let allPtsBeyond x = all (\y -> 1-epsilon <= (x % y)) adjPts in
  -- | From the candidate hyperplanes, choose the ones that truly separate the
  -- classes.
  let corners = filter allPtsBeyond candidates in 
  if null corners then error "No corners for analysis." else
    -- Now from those, choose the one with the tightest fit.
    let (b : w, minWSqr) = minimize (sqr . tail) corners in
    normalize (LC (w, b))

normalize :: Floating a => LC a -> LC a
--normalize = normalizeByVector
normalize = normalizeByCutoff

normalizeByVector :: Floating a => LC a -> LC a
normalizeByVector (LC (w, b)) = LC ((1/mag) *! w, b/mag)
    where mag = size w

normalizeByCutoff :: Fractional a => LC a -> LC a
normalizeByCutoff (LC (w, b)) =
    LC ((1/b) *! w, b/b)


type Dataset = [(Bool, [Double])]
example :: Dataset
example = [(False, [0,0,1]), (True, [1,0,0])]
example2 :: Dataset
example2 = [(False, [0,0,1]), (True, [1,0,0]), (False, [0,1,1]), (True, [0,0,-1])]
example3 :: Dataset
example3 = [(True, [2,1]), (False, [3,2]), (True, [6,2]),
            (False, [5,5]), (True, [6,3])]

instance (Random a, Random b, Random c)
    => Random (a, b, c) where
    randomR ((lo1, lo2, lo3), (hi1, hi2, hi3)) gen = 
        let (x, gen') = randomR (lo1, hi1) gen in
        let (y, gen'') = randomR (lo2, hi2) gen' in
        let (z, gen''') = randomR (lo3, hi3) gen'' in
        ((x, y, z), gen''')
    random gen = 
        let (x, gen') = random gen in
        let (y, gen'') = random gen' in
        let (z, gen''') = random gen'' in
        ((x, y, z), gen''')

-- | Randomly generated dataset that has nothing in the margin z \in
-- (0.4999, 0.5001) and tags everything based on the sign of z. Note:
-- Doesn't work! There is a bug somewhere!
genExample :: Int -> IO Dataset
genExample n = replicateM n (do (x,y,z) <- randomIO
                                let tag = z < 0
                                let zoffset = if tag then -0.001 else 0.001
                                return (tag, [x,y,z+zoffset]))

-- | Given a dataset, find a classifier, then run the classifier on
-- all the points and see how many classifications are (good, bad).
test_example :: Dataset -> (Int, Int)
test_example dataset =
  let clfier = margin dataset in
  let results = map (\(tag, vec) -> (tag, classify clfier vec)) dataset in
  let (good, bad) = partition (\(x,y) -> x==y) results in
  (length good, length bad)

---- List utilities: --

-- (TODO: unused).
unzipWith :: (a -> (s, t)) -> [a] -> ([s], [t])
unzipWith f xs = unzip $ map f xs

-- | All the values @(a, b)! where @a@ and @b@ come from the given list,
-- and @a@ precedes @b@ in that list.
-- (TODO: unused)
allPairs :: [a] -> [(a, a)]
allPairs [] = []
allPairs (x:ys) = map (\y -> (x,y)) ys ++ allPairs ys

-- | Generate all subsets of a given size from a given list.
-- Running time: T(n,xs) = n + T(n-1,xs) + T(n, xs-1).
-- Produces (n C xs) = xs! / (xs-n)!n! outputs.
combos :: Int -> [a] -> [[a]]
combos 0 _ = [[]]
combos n xs | length xs < n  =  []
combos n (x:xs) = map (x :) (combos (n-1) xs) ++ combos n xs

joinMaybes :: MonadPlus m => m (Maybe a) -> m a
joinMaybes = (>>= maybe mzero return)

-- | Return the "graph" of a function @f@, restricted to inputs @xs@.
-- That is, all the pairs @(x, f x)@ for inputs in @xs@.
graph :: (a -> b) -> [a] -> [(a, b)]
graph f xs = [(x, f x) | x <- xs]

minimize :: Ord b => (a -> b) -> [a] -> (a, b)
minimize f [] = error "Can't minimize over no data points"
minimize f xs = foldr1 (\(x1, y1) (x2, y2) ->
                        if y1 > y2 then (x1, y1) else (x2, y2))
                (graph f xs)

-- | Return all the rotations of @xs@. (All values of @roll n xs@ for
-- n <- [0..length n-1].)]
rollings :: [a] -> [[a]]
rollings xs = take n $ map (take n) $ tails $ cycle xs
    where n = (length xs)
