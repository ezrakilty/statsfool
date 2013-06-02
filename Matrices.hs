module Matrices where

import Control.Monad (MonadPlus, mzero, msum)
import Test.HUnit

import Debug.Trace

---- Vector operations: --

-- | The type of a vector over field @a@.
type Vector a = [a]

(%) :: Num a => Vector a -> Vector a -> a
u % v = sum $ zipWith (*) u v  -- ^ dot product

(*!) :: Num a => a -> Vector a -> Vector a
infixr 9 *!
k *! u = map (*k) u  -- ^ scalar multiplication

(+!) :: Num a => Vector a -> Vector a -> Vector a
infixr 8 +!
u +! v = zipWith (+) u v  -- ^ vector addition

(-!) :: Num a => Vector a -> Vector a -> Vector a
infixr 8 -!
u -! v = zipWith (-) u v  -- ^ vector subtraction

-- | The negation of a vector.
neg :: Num a => Vector a -> Vector a
neg = map negate

-- | The dimension of a vector, or number of its elements.
dimension :: [a] -> Int
dimension x = length x

-- | The magnitude of a vector.
size :: Floating a => Vector a -> a
size x = sqrt (x % x)

-- | The squared-magnitude of a vector (self-dot-product).
sqr :: Num a => Vector a -> a
sqr x = x % x

-- | Is the given vector a unit vector? Allows a tolerance of episilon
-- of the squared magnitude.
isUnitVec :: Vector Double -> Bool
isUnitVec vec = abs(sqr vec - 1) <= epsilon

-- | Is the given vector exactly zero?
isZeroVec :: (Eq a, Num a) => Vector a -> Bool
isZeroVec = all (==0)

-- | Is the given vector either a unit-vector or a zero-vector? This
-- is exact; it allows for no tolerance.
isZeroOrUnitBasisVec :: (Eq a, Num a) => Vector a -> Bool
isZeroOrUnitBasisVec [] = True
isZeroOrUnitBasisVec (0:xs) = isZeroOrUnitBasisVec xs
isZeroOrUnitBasisVec (1:xs) = isZeroVec xs
isZeroOrUnitBasisVec (_: _) = False

---- Matrix operations: --

-- | The type of a matrix over field @a@.
type Matrix a = [[a]]

-- | The identity matrix, size 3x3.
id3 :: Matrix Double
id3 = [[1,0,0], [0,1,0], [0,0,1]]

-- | Transpose a matrix.
transpose :: Matrix a -> Matrix a
transpose [] = [[]]
transpose ([]:_) = []
transpose m = map head m : transpose (map tail m)

-- | Return a list of the columns of a matrix
cols :: Matrix a -> [Vector a]
cols = transpose

-- | Return a list of the rows of a matrix
rows :: Matrix a -> [Vector a]
rows x = x

-- | The number of rows (height) of a matrix.
nRows :: Matrix a -> Int
nRows = length . rows
-- | The number of columns (width) of a matrix.
nCols :: Matrix a -> Int
nCols = length . cols
-- | The number of rows (height) of a matrix.
height :: Matrix a -> Int
height = nRows
-- | The number of columns (width) of a matrix.
width :: Matrix a -> Int
width = nCols

-- | Matrix multiplication, naively implemented.
matMult :: Num a => Matrix a -> Matrix a -> Matrix a
matMult m n | nCols m == nRows n
  = [[mRow % nCol | nCol <- cols n] | mRow <- rows m]
matMult m n = error "Row/column size mismatch in matrix multiplication."

-- | Infix operator for matrix multiplication.
(%%) :: Num a => Matrix a -> Matrix a -> Matrix a
(%%) = matMult

-- | @transform m v@ applies the linear transformation designated by
-- matrix @m@ to the vector @v@.
transform :: Num a => Matrix a  -> Vector a -> Vector a
transform m v = fromColVec (m %% (toColVec v))

-- | Given two matrices @l@ and @r@, append them horizontally.
paste :: Matrix a -> Matrix a -> Matrix a
paste = zipWith (++)

---- Numerical mumbojumbo: --

epsilon :: Double
epsilon = 0.0000001 -- ^ tolerance for round-off errors

-- | Are the given vectors within @episilon@ of each other? (By
-- squared magnitude.)
closeTo :: Vector Double -> Vector Double -> Bool
closeTo u v = sqr (u -! v) <= epsilon

-- | @annhilate u v@ subtracts from @v@ the multiple of @u@ that produces a
-- zero in the first coordinate of @u@
annhilate :: Fractional a => Vector a -> Vector a -> Vector a
annhilate (u@(x:_)) (v@(y:_)) =
    let x:_ = u ; y:_ = v in
      v -! ((y / x) *! u)

-- | Given two vectors @u@, @v@, representing rows in a Gaussian elimination
-- matrix, what factor should we apply to @u@ to make it
-- eliminate the first non-zero component of @v@.
gaussFrac :: (MonadPlus m, Eq a, Fractional a) =>
             [a] -> [a] -> m a
gaussFrac u v = msum (zipWith f u v) 
  where f 0 _ = mzero
        f a b = return (b / a)

-- | Return the first element of a list which is not 0.
firstNonZero :: (Eq a, Num a) => [a] -> a
firstNonZero [] = error "no non-zero element in vector passed to firstNonZero"
firstNonZero (0:xs) = firstNonZero xs
firstNonZero (x:_) = x

-- | Normalize a vector for Gaussian elimination: scale it so that its
-- first component is approximately 1. (Only works on vectors with
-- nonzero first component.)
gaussNormalize :: (Eq a, Fractional a) => [a] -> [a]
gaussNormalize row | head row /= 0 = (1/head row) *! row

-- | Convert a matrix to upper triangular form, with ones in the
-- leading non-zero positions, by Gaussian elimination.
upperTriang :: (Eq a, Floating a) => Matrix a -> Matrix a
upperTriang [] = []
upperTriang (matrix@([]:_)) = matrix
upperTriang (row@(0:_):rows) = row : upperTriang rows
upperTriang (row:rows) = let row' = gaussNormalize row
                             rows' = map (annhilate row) rows in
                         row' : (map (0:) $ upperTriang (map tail rows'))

-- | Apply a step of Gaussian elimination to two vectors: Given @u@
-- and @v@, subtract a multiple of @u@ from @v@ so that @v@ has no
-- component in the direction of @u@'s first nonzero component.
gaussStep :: (Eq a, Floating a) => [a] -> [a] -> [a]
gaussStep u v = 
    case k of Nothing -> v
              Just  k -> v -! k*!u
    where k = gaussFrac u v

-- | @upperToEchelon@ takes an upper-triangular matrix and returns an
-- equivalent echelon matrix. ("Equivalent" as a system of linear
-- equations.)
upperToEchelon :: (Eq a, Floating a) => Matrix a -> Matrix a
upperToEchelon [] = []
upperToEchelon [row] = [row]
upperToEchelon (row:rows) = row : upperToEchelon (map (gaussStep row) rows)

-- | Is the given matrix an echelon matrix? An echelon matrix is one
-- where each row has a leading one, preceded by all zeroes, and the
-- other entries in that column are all zero. Furthermore, the 1 in
-- the i-th row should be in a column less than the column of the 1 in
-- the j-th row, for all i and j.
isEchelon' :: Matrix Double -> Bool
isEchelon' m = and $ take n $ [isZeroOrUnitBasisVec c | c <- cols m]
    where n = height m

-- | Is the given matrix an echelon matrix? (Alternate implementation.)
isEchelon :: Matrix Double -> Bool
isEchelon m = all (isUnitVec) $ take (nRows m) $ 
               filter (not . isZeroVec) (cols m)

-- | Given two lists of values, where the first one contains a 1,
-- return the element of the second list in the position corresponding
-- to that 1; throw error if there is no such 1.
genericSelect :: (Eq a, Num a) => [a] -> [b] -> b
genericSelect [] _ = error "No 1 in first argument to select"
genericSelect (1:xs) (y:ys) = y
genericSelect (_:xs) (_:ys) = genericSelect xs ys

select :: [Double] -> [b] -> b
select [] _ = error "No 1 in first argument to select"
select (x:xs) (y:ys) | abs(x-1)<=epsilon = y
select (_:xs) (_:ys) = select xs ys

-- | Like @select@, but works with a list of floating point numbers,
-- using the first position in the first list that is within @epsilon@ of 1.
select1 :: Double -> Double -> Maybe Double
select1 x y | abs(x-1) < epsilon = Just y
            | otherwise = Nothing

-- | Given an upper-echelon representation of a system of linear
-- equations with a unique solution @v@, return @Just v@; if there
-- is not a unique solution, return @Nothing@.
extractSol :: Matrix Double -> Maybe (Vector Double)
extractSol m = sequence $
    for coeffs $ \u ->
        if isUnitVec u then
             return (u `select` vals) 
        else fail "non-unit column ..."
    where m' = transpose m
          vals = last m'
          coeffs = init m'

-- | Using Gaussian elimination, solve a system of linear equations.
-- Returns the solution set represented as an upper-echelon matrix.
gaussianElim :: Matrix Double -> Matrix Double
gaussianElim m =
    reverse $ upperToEchelon $ reverse $ upperTriang m

-- | Return @Just v@ if @v@ is the unique solution vector for @m@. If there
-- is no unique solution, return @Nothing@.
gaussianElimUniq :: Matrix Double -> Maybe (Vector Double)
gaussianElimUniq m
    | width m <= height m + 1  =  extractSol (gaussianElim m)
    | otherwise = fail "Linear system is under-determined."

-- | Turn a vertical 1xn matrix into a row vector.
fromColVec :: Matrix a -> Vector a
fromColVec = map head

-- | Turn a row vector into a vertical 1xn matrix.
toColVec :: Vector a -> Matrix a
toColVec = map return

-- | Check that @soln@ is a solution to the linear system specified by
-- (@m@, @v@). (Within @epsilon@.)
solves :: Matrix Double -> Vector Double -> Vector Double -> Bool
solves m v soln = 
    closeTo v (fromColVec (matMult m (toColVec soln)))
--    closeTo v $ fromColVec $ matMult m $ toColVec soln
--    closeTo v . fromColVec . matMult m . toColVec $ soln

-- | @m <|> v@ pastes the vector @v@, in column form, onto the right
-- side of @m@. Alternatively: Given a matrix @m@ of variable
-- coefficients and a vector @v@ of constants to which we equate
-- them, return a linear system (as a matrix)--ready to send to
-- @gaussianElimUniq@.
(<|>) :: Matrix a -> Vector a -> Matrix a
m <|> v = paste m (toColVec v)

---- Tests: --

-- | Test (with @solves@) that @guassianElimUniq@ produces a solution
-- to the system (@m@, @v@).
prop_gaussianElim :: Matrix Double -> Vector Double -> Bool
prop_gaussianElim m v = case (gaussianElimUniq (m <|> v)) of
                        Nothing -> False
                        Just soln -> solves m v soln

test_gaussianElim :: Test
test_gaussianElim = test [assert $ prop_gaussianElim [[2, 4], [3,-1]] [6,5],
                          assert $ prop_gaussianElim [] []]

---- More matrix operations: --

-- | The list of all NW-SE diagonals of a matrix.
seDiags :: Matrix a -> [Vector a]
seDiags m = transpose $ zipWith roll [0..] m

-- | The list of all NE-SW diagonals of a matrix.
neDiags :: Matrix a -> [Vector a]
neDiags m = seDiags $ reverse m

-- | The determinant of a matrix.
determinant :: Num a => Matrix a -> a
determinant m = sum (map product (seDiags m))
                - sum (map product (neDiags m))

-- | The determinant of a matrix.
det :: Num a => Matrix a -> a
det = determinant

---- List utilities --
-- | Like @map@ but affords a nicer syntax, putting the body of the
-- function last.
for :: [a] -> (a->b) -> [b]
for xs f = map f xs

-- | Rotate the members of @xs@ (treated as a circular list) by @n@ places.
roll :: Int -> [a] -> [a]
roll n xs = take len $ drop n $ cycle xs
    where len = length xs
