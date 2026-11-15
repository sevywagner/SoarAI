import { createBrowserRouter } from "react-router";
import Header from "./layout/Header";
import ProcessMassValue from "../pages/ProcessMassValue";
import Graphs from "../pages/Graphs";
import About from "../pages/About";
import Error from "../pages/Error";

const router = createBrowserRouter([
    {
        path: '/mass-spec-web-ui',
        element: <Header />,
        errorElement: <Error />,
        children: [
            {
                path: '/mass-spec-web-ui/',
                element: <ProcessMassValue />
            },
            {
                path: '/mass-spec-web-ui/graphs',
                element: <Graphs />
            }
        ]
    }
]);

export default router;
