import java.io.File

import org.apache.commons.io.FileUtils

import scala.sys.process._

object main extends App {
  val testsDir: String = "../exported_tests"
  val matrixmul: String = "../cmake-build-debug/matrixmul"
  val dir = new File(testsDir)
  private val files: Array[File] = dir.listFiles.filter(_.isFile)

  files.filter(_.getName.startsWith("result")).sortWith((f1, f2) => f1.getName < f2.getName).foreach(result => {
    val name = """result_(\d+)_000(\d+)_(\d+)_(\d+)""".r
    result.getName match {
      case name(x, y, z, a) => {
        val sparseMatrix = files.filter(_.getName.equals(s"sparse05_000${y}_$z")).head
        val denseMatrix = files.filter(_.getName.equals(s"matrix01_000${y}_$a")).head

        val myResultFile = new File(s"../cmake-build-debug/outputs_tests/${result.getName}")
        myResultFile.createNewFile()

        for(i <- 1 to 2) {
          val runCommand = s"mpiexec -n $i $matrixmul -f ${sparseMatrix.getCanonicalPath} -s $a -c 1 -e $x"
          val ret = runCommand #> myResultFile !

          println(s"Test ${result.getName} \n$runCommand")
          val failure: String = s" Failure :( "
          val success: String = s" Success "
          if (FileUtils.contentEquals(myResultFile, result)) {
            println(success)
          } else {
            println(failure)
            //          "cat" #< myResultFile !
          }
          println()
        }


      }
    }
  })
}
