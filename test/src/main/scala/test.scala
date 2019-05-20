import java.io.File
import java.net.URL

import scala.sys.process._
import scala.collection.JavaConverters._

object test extends App{

  private val strings: Stream[String] = "cat" #< new File("file") lineStream;

  strings.

  print(("cat" #< strings) !! );


}
