import java.io.File

import geny.Generator.from

import scala.math.abs
import scala.util.Random
import scala.collection.immutable.List._
import scala.sys.process._

object generateSparse extends App {

  def printToFile(f: java.io.File)(op: java.io.PrintWriter => Unit) {
    val p = new java.io.PrintWriter(f)
    try { op(p) } finally { p.close() }
  }

  val n = 24*24*12
//  val n = 10
  val maxNonzero=42
  val nonZeros = 1000000
//  val nonZeros = 10
  val (valuesBegin, valuesEnd) = (0,1.0)

  var random = new Random()
  println("before firstLine")
  val firstLine = (1 to nonZeros).map(_ => random.nextFloat() * (valuesEnd - valuesBegin) + valuesBegin)
    .map(a => "%.5f".format(a)).foldLeft(new StringBuilder())((b,s) => b.append(s).append(" ")).toString()
  println("after firstLine")
  private val extentsVector = (nonZeros :: 0 :: List.fill(n-1)(abs(random.nextInt()) % (nonZeros + 1))).sorted.to[Vector]
  val extents = extentsVector.reduce[Any]((a, b) => s"${a toString} ${b toString}")
  println("after extents")
  val indices = (0 to n - 1).zip(1 to n)
    .flatMap{ case (a,b) => {
      val range = (extentsVector(a) until extentsVector(b)).toList
      random.shuffle((0 until n).toList).take(range.length)
    }}.map(a => a.toString).foldLeft(new StringBuilder())((b,s) => b.append(s).append(" ")).toString()
  println("after indices")

  private val testFile = new File("../tests/2")
  testFile.createNewFile()
  printToFile(testFile) { p =>
    p.println(s"$n $n $nonZeros $maxNonzero")
    p.println(firstLine)
    p.println(extents)
    p.println(indices)
  }

}
