import Control.Monad
import Control.Monad.ST
import Data.Array.ST
import Data.STRef
import Numeric.LinearAlgebra 
import Data.STRef
import Data.Time.Clock
import Data.Maybe

gaussianElimination :: Matrix Double -> (Matrix Double, Double)
gaussianElimination mat = runST $ do
  let n = rows mat
  let elements = concat $ toLists mat
  stMatrix <- newListArray ((0, 0), (n-1, n-1)) elements :: ST s (STArray s (Int, Int) Double)
  pivotRef <- newSTRef (-1)
  signRef <- newSTRef 1

  forM_ [0..n-1] $ \col -> do
    -- pivot
    forM_ [col..n-1] $ \row -> do
      val <- readArray stMatrix (row, col)
      when (val /= 0) $ do
        writeSTRef pivotRef row
        return ()

    pivot <- readSTRef pivotRef

    if pivot == -1
      then writeSTRef signRef 0
      else do
        -- swap
        when (pivot /= col) $ do
          modifySTRef signRef negate
          forM_ [0..n-1] $ \i -> do
            temp1 <- readArray stMatrix (col, i)
            temp2 <- readArray stMatrix (pivot, i)
            writeArray stMatrix (col, i) temp2
            writeArray stMatrix (pivot, i) temp1

        -- elimination
        pivotElem <- readArray stMatrix (col, col)
        forM_ [col+1..n-1] $ \row -> do
          rowElem <- readArray stMatrix (row, col)
          forM_ [col..n-1] $ \i -> do
            val <- readArray stMatrix (row, i)
            pivotVal <- readArray stMatrix (col, i)
            writeArray stMatrix (row, i) (val - rowElem * pivotVal / pivotElem)

  sign <- readSTRef signRef
  elements <- getElems stMatrix
  return (reshape n (fromList elements),sign)


determinant :: Matrix Double -> Double
determinant mat =
  let (filteredMat, sign) = gaussianElimination mat
      dett = product $ toList $ takeDiag filteredMat
  in dett * sign

calcRank :: Matrix Double -> Int
calcRank mat =
  let (filteredMat, _) = gaussianElimination mat
      nonZeroRows = filter (not . all (\x -> abs x < 1e-10)) $ toLists filteredMat
  in length nonZeroRows  


inverseMatrix :: Matrix Double -> Maybe (Matrix Double)
inverseMatrix mat = runST $ do 
  let n = rows mat
  let c = 2*n
  let identMatrix = ident n
  let mergedMatrix = fromBlocks [[mat,identMatrix]]
  
  let elements = concat $ toLists mergedMatrix
  stMatrix <- newListArray ((0, 0), (n-1, c-1)) elements :: ST s (STArray s (Int, Int) Double)
  pivotRef <- newSTRef (-1)

  isSingular <- newSTRef False  -- Flag to detect singular matrix

  forM_ [0..n-1] $ \col -> do
    forM_ [col..n-1] $ \row -> do
      val <- readArray stMatrix (row, col)
      when (val /= 0) $ writeSTRef pivotRef row

    pivot <- readSTRef pivotRef

    if pivot == -1
      then writeSTRef isSingular True
      else do
        when (pivot /= col) $ do
          forM_ [0..c-1] $ \i -> do
            temp1 <- readArray stMatrix (col, i)
            temp2 <- readArray stMatrix (pivot, i)
            writeArray stMatrix (col, i) temp2
            writeArray stMatrix (pivot, i) temp1
        
        factor <- readArray stMatrix (col, col)
        forM_ [col..c-1] $ \i -> do
          val <- readArray stMatrix (col,i)
          let normalizeVal = val / factor
          writeArray stMatrix (col,i) normalizeVal

        forM_ [0..n-1] $ \row -> when (row /= col) $ do
          rowElem <- readArray stMatrix (row, col)
          forM_ [col..c-1] $ \i -> do
            val <- readArray stMatrix (row, i)
            pivotVal <- readArray stMatrix (col, i)
            writeArray stMatrix (row, i) (val - rowElem * pivotVal)

  singular <- readSTRef isSingular
  if singular
    then return Nothing
    else do
      invElements <- forM [0..n-1] $ \row -> 
        forM [n..c-1] $ \col -> readArray stMatrix (row, col)

      -- Convert to Vector and reshape
      let flatInvElements = concat invElements
      return $ Just (reshape n (fromList flatInvElements))


-- Function to create a matrix from a text file
buildMatrixFromFile :: FilePath -> IO (Matrix Double)
buildMatrixFromFile filePath = do
    content <- readFile filePath
    let rows = map (map read . words) (lines content)
    return $ fromLists rows

printMatrixElements :: Matrix Double -> IO ()
printMatrixElements mat = do
    let rows = toLists mat            
    let firstRow = take 5 (head rows)  
    let lastRow = take 5 (last rows)   
    putStrLn "First 5 elements of the first row:"
    print firstRow
    putStrLn "First 5 elements of the last row:"
    print lastRow

-- main :: IO ()
-- main = do
--   tt_start <- getCurrentTime
--   putStrLn "Inverse Matrix : "
--   matrix <- buildMatrixFromFile "m1.txt"
--   case inverseMatrix matrix of
--     Just invMatrix -> printMatrixElements invMatrix
--     Nothing -> putStrLn "Matrix is not invertible"
--   tt_end <- getCurrentTime
--   putStrLn "Time :"
--   print (diffUTCTime tt_end tt_start)

-- main :: IO ()
-- main = do 
--     mat <- buildMatrixFromFile "Rank.txt"
--     tt_start <- getCurrentTime
--     putStrLn "Rank :  "
--     let rank = calcRank mat 
--     print rank
--     tt_end <- getCurrentTime
--     putStrLn "Time"
--     print (diffUTCTime tt_end tt_start)

-- main :: IO ()
-- main = do 
--     mat <- buildMatrixFromFile "Determinant.txt"
--     tt_start <- getCurrentTime
--     putStrLn "Rank :  "
--     let rank = determinant mat 
--     print rank
--     tt_end <- getCurrentTime
--     putStrLn "Time"
--     print (diffUTCTime tt_end tt_start)