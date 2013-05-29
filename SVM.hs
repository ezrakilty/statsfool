-- | A specification for a Support Vector Machine algorithm. VERY SLOW!

-- Copyright Ezra Cooper, circa 2008.

{-
  TODO:
    * Choose the candidates more cleverly.
    * Use Data.Vector for matrix representations.  
    * Examine the Gaussian algorithm. Faster versions?
 -}

module SVM where

import Test.HUnit
import Data.List (tails)
import Control.Monad (MonadPlus, mplus, mzero, msum)

---- Vector operations: --

type Vector a = [a]

u % v = sum $ zipWith (*) u v  -- ^ dot product

infixr 9 *!
k *! u = map (*k) u  -- ^ scalar multiplication

infixr 8 +!
u +! v = zipWith (+) u v  -- ^ vector addition

infixr 8 -!
u -! v = zipWith (-) u v  -- ^ vector subtraction

neg = map negate

dimension x = length x

size x = sqrt (x % x)

sqr x = x % x

isUnitVec vec = abs(sqr vec - 1) <= epsilon

isZeroVec = all (==0)

isZeroOrUnitBasisVec [] = True
isZeroOrUnitBasisVec (0:xs) = isZeroOrUnitBasisVec xs
isZeroOrUnitBasisVec (1:xs) = isZeroVec xs
isZeroOrUnitBasisVec (_: _) = False

---- Matrix operations: --

type Matrix a = [[a]]

id3 :: Matrix Double
id3 = [[1,0,0], [0,1,0], [0,0,1]]

transpose [] = [[]]
transpose ([]:_) = []
transpose m = map head m : transpose (map tail m)

cols = transpose
rows x = x

nRows = length . rows
nCols = length . cols
height = nRows
width = nCols

matMult m n | nCols m == nRows n
  = [[mRow % nCol | nCol <- cols n] | mRow <- rows m]
matMult m n = error "row/column mismatch in matrix mult, e.g, non-conformable matrices"

(%%) = matMult

-- | @transform m v@ applies the linear transformation designated by
-- matrix @m@ to the vector @v@.
transform m v = m %% (toColVec v)

---- Numerical mumbojumbo: --

epsilon :: Double
epsilon = 0.0000001 -- tolerance for round-off errors

closeTo u v = sqr (u -! v) <= epsilon

-- | @annhilate u v@ subtracts from @v@ the multiple of @u@ that produces a
-- zero in the first coordinate of @u@
annhilate (u@(x:_)) (v@(y:_)) =
    let x:_ = u ; y:_ = v in
      v -! ((y / x) *! u)

gaussFrac u v = msum (zipWith f u v) 
  where f 0 _ = mzero
        f a b = return (b / a)

-- | Return the first element of a list which is not 0.
firstNonZero [] = error "no non-zero element in vector passed to firstNonZero"
firstNonZero (0:xs) = firstNonZero xs
firstNonZero (x:_) = x

gaussNormalize row | head row /= 0 = (1/head row) *! row

-- | Convert a matrix to upper triangular form, with ones in the
-- leading non-zero positions, by Gaussian elimination
upperTriang :: (Eq a, Floating a) => Matrix a -> Matrix a
upperTriang [] = []
upperTriang (matrix@([]:_)) = matrix
upperTriang (row@(0:_):rows) = row : upperTriang rows
upperTriang (row:rows) = let row' = gaussNormalize row
                             rows' = map (annhilate row) rows in
                         row' : (map (0:) $ upperTriang (map tail rows'))

gaussStep :: (Eq a, Floating a) => [a] -> [a] -> [a]
gaussStep u v = 
    case k of Nothing -> v
              Just k -> v -! k*!u
    where k = gaussFrac u v

-- | @upperToEchelon@ takes an upper-triangular matrix and returns an
-- equivalent echelon matrix.
upperToEchelon :: (Eq a, Floating a) => Matrix a -> Matrix a
upperToEchelon [] = []
upperToEchelon [row] = [row]
upperToEchelon (row:rows) = row : upperToEchelon (map (gaussStep row) rows)

-- | An echelon matrix is one where each row has a leading one, preceded
-- by all zeroes, and the other entries in that column are all zero.
-- Furthermore, the 1 in the i-th row should be in a column less than
-- the column of the 1 in the j-th row, for all i and j.
isEchelon' m = and $ take n $ [isZeroOrUnitBasisVec c | c <- cols m]
    where n = height m

isEchelon m = all (isUnitVec) $ take (nRows m) $ 
               filter (not . isZeroVec) (cols m)

-- | Given two lists of values, where the first one contains a 1,
-- return the element of the second list in the position corresponding
-- to that 1; throw error if there is no such 1.
select [] _ = error "No 1 in first argument to select"
select (1:xs) (y:ys) = y
select (_:xs) (_:ys) = select xs ys

-- | Like @select@, but works with a list of floating point numbers,
-- using the first position in the first list that is within @epsilon@ of 1.
select1 x y | abs(x-1) < epsilon = Just y
            | otherwise = Nothing

-- | Given an upper-echelon representation of a system of linear
-- equations with a unique solution v, return Just v; if there
-- is not a unique solution, return Nothing.
extractSol :: Matrix Double -> Maybe (Vector Double)
extractSol m = sequence $
    for coeffs $ \u ->
        if isUnitVec u then
             return (u `select` vals) 
        else fail "non-unit column ..."
    where m' = transpose m
          vals = last m'
          coeffs = init m'

-- | Uses gaussian elimination to solve a system of linear equations.
-- Returns the solution set represented as an upper-echelon matrix.
gaussianElim :: Matrix Double -> Matrix Double
gaussianElim m =
    reverse $ upperToEchelon $ reverse $ upperTriang m

-- | Return @Just v@ if @v@ is the unique solution vector for @m@. If there
-- is no unique solution, return @Nothing@.
gaussianElimUniq :: Matrix Double -> Maybe (Vector Double)
gaussianElimUniq m | width m <= height m + 1  =  extractSol (gaussianElim m)
                   | otherwise = fail "Matrix is not shaped for a unique solution."

-- | Turn a vertical 1xn matrix into a row vector.
fromColVec = map head

-- | Turn a row vector into a vertical 1xn matrix.
toColVec = map return :: [a] -> [[a]]

solves m v soln = 
    closeTo v (fromColVec $ matMult m (toColVec soln))

makeSystem weights vec = zipWith snoc weights vec

-- | @m <|> v@ pastes the vector @v@, in column form, onto the right side of @m@.
m <|> v = paste m (toColVec v)

---- Tests: --

prop_gaussianElim m v = case (gaussianElimUniq (m <|> v)) of
                        Nothing -> False
                        Just soln -> solves m v soln

test_gaussianElim = test [assert $ prop_gaussianElim [[2, 4], [3,-1]] [6,5],
                          assert $ prop_gaussianElim [] []]

---- More matrix operations: --

-- | The list of all NW-SE diagonals of a matrix.
seDiags m = transpose $ zipWith roll [0..] m

-- | The list of all NE-SW diagonals of a matrix.
neDiags m = seDiags $ reverse m

-- | The determinant of a matrix.
determinant :: Num a => Matrix a -> a
determinant m = sum (map product (seDiags m))
                - sum (map product (neDiags m))

-- | The determinant of a matrix.
det = determinant

-- linear classifier
data LC a  =  LC(Vector a, a)
    deriving Show

classify :: (Ord a, Floating a) => LC a -> Vector a -> Bool
classify (LC(w, b)) x = w % x - b <= 0

-- | Find the maximal-margin hyperplane
--margin :: [(Bool, Vector a)] -> Vector a
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
    LC (w, b)

type Dataset = [(Bool, [Double])]
example :: Dataset
example = [(False, [0,0,1]), (True, [1,0,0])]
example2 :: Dataset
example2 = [(True, [2,1]), (False, [3,2]), (True, [6,2]),
            (False, [5,5]), (True, [6,3])]

---- List utilities: --

unzipWith f xs = unzip $ map f xs

allPairs [] = []
allPairs (x:ys) = map (\y -> (x,y)) ys ++ allPairs ys

-- | Generate all subsets of a given size from a given list.
-- Running time: T(n,xs) = n + T(n-1,xs) + T(n, xs-1).
-- Produces (n C xs) = xs! / (xs-n)!n! outputs.
combos 0 _ = [[]]
combos n xs | length xs < n  =  []
combos n (x:xs) = map (x :) (combos (n-1) xs) ++ combos n xs

joinMaybes :: MonadPlus m => m (Maybe a) -> m a
joinMaybes = (>>= maybe mzero return)

graph f xs = [(x, f x) | x <- xs]

minimize :: Ord b => (a -> b) -> [a] -> (a, b)
minimize f [] = error "Can't minimize over no data points"
minimize f xs = foldr1 (\(x1, y1) (x2, y2) ->
                        if y1 > y2 then (x1, y1) else (x2, y2))
                (graph f xs)

paste = zipWith (++)

snoc xs x = xs ++ [x]

rollings xs = take n $ map (take n) $ tails $ cycle xs
    where n = (length xs)

roll n xs = take len $ drop n $ cycle xs
    where len = length xs

-- | Like @map@ but affords a nicer syntax, putting the body of the
-- function last.
for xs f = map f xs
