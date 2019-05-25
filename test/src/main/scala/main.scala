import java.io.File
import java.util.Scanner
import scala.math

import org.apache.commons.io.FileUtils

import scala.io.Source
import scala.sys.process._

object main extends App {
  def filesEqual(file1 : File, file2 : File) : Boolean= {
    var result : Boolean = true
    var scanner1 = new Scanner(file1)
    var scanner2 = new Scanner(file2)
    while(scanner1.hasNextDouble) {
      if (!scanner2.hasNextDouble)
        result = false
      val d1: Double = scanner1.nextDouble()
      val d2: Double = scanner2.nextDouble()
      val epsilon : Double = 0.000012
      if (scala.math.abs(d1 - d2) > epsilon)
        result = false
    }
    if(scanner2.hasNextDouble) {
      result = false
    }

    return result
  }
  val testsDir: String = "../exported_tests"
  val matrixmul: String = "../cmake-build-debug/matrixmul"
  val dir = new File(testsDir)
  private val files: Array[File] = dir.listFiles.filter(_.isFile)

  var successes : Int = 0
  var failures : Int = 0

  files.filter(_.getName.startsWith("result")).sortWith((f1, f2) => f1.getName < f2.getName).foreach(result => {
    val name = """result_(\d+)_000(\d+)_(\d+)_(\d+)""".r
    result.getName match {
      case name(x, y, z, a) => {
        val sparseMatrix = files.filter(_.getName.equals(s"sparse05_000${y}_$z")).head
        val denseMatrix = files.filter(_.getName.equals(s"matrix01_000${y}_$a")).head

        val myResultFile = new File(s"../cmake-build-debug/outputs_tests/${result.getName}")
        myResultFile.createNewFile()
        val paramteresCombinations = List((2, 2), (4, 2), (6, 2), (6, 3))
        List((2,1)).foreach {
          case (n, c) => {
            val runCommand = s"mpiexec -n $n $matrixmul -f ${sparseMatrix.getCanonicalPath} -s $a -c $c -e $x -v -i"
            val ret = runCommand #> myResultFile !

            println(s"Test ${result.getName} \n$runCommand")
            val failure: String = s" Failure :( "
            val success: String = s" Success "
            if (filesEqual(myResultFile, result)) {
              println(success)
              successes += 1
            } else {
              println(failure)
              failures += 1
              //          "cat" #< myResultFile !
            }
            println()
          }
        }


      }
    }
  })
  println(s"Failures $failures Successes $successes")
}
