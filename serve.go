package coverageplot

import (
	"github.com/gorilla/mux"
	"github.com/meatballhat/negroni-logrus"
	"github.com/urfave/negroni"
	"log"
	"net/http"
)

func Serve() *http.Server {

	srv := &http.Server{Addr: ":8080"}
	r := mux.NewRouter()

	//put your regular handlers here

	//then comes root handler
	// r.HandleFunc("/", homePageHandler)
	//if a path not found until now, e.g. "/image/tiny.png"
	//this will look at "./public/image/tiny.png" at filesystem
	r.PathPrefix("/").Handler(http.FileServer(http.Dir("../public/")))

	// Set up a request logger, useful for debugging
	n := negroni.New()
	n.Use(negronilogrus.NewMiddleware())
	n.UseHandler(r)

	http.Handle("/", n)

	go func() {
		if err := srv.ListenAndServe(); err != nil {
			// cannot panic, because this probably is an intentional close
			log.Printf("Httpserver: ListenAndServe() error: %s", err)
		}
	}()

	return srv
}
