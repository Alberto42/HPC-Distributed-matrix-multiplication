import java.io.File
import java.net.URL

import scala.sys.process._
import scala.collection.JavaConverters._

object test extends App{


  List((1,2)).foreach(x => {
    x match {
      case (a,b) => println(s"dupa $a dupa $b")
    }
  })


}
